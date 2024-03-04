# ------------------------------------------------------------------------------
# 
# 01_prepare_monthly_data
# 
# Read ISIMIP NetCDF rasters and transform them to R data frames describing
# monthly values (columns) in grid cells (rows).
# 
# Analysis script for publication:
# 
# "Notable shifts beyond pre-industrial streamflow and soil moisture conditions
# transgress the planetary boundary for freshwater change"
#
# Miina Porkka, Vili Virkki, Lan Wang-Erlandsson, Dieter Gerten, Tom Gleeson,
# Chinchu Mohan, Ingo Fetzer, Fernando Jaramillo, Arie Staal, Sofie te Wierik,
# Arne Tobian, Ruud van der Ent, Petra Döll, Martina Flörke, Simon N. Gosling,
# Naota Hanasaki, Yusuke Satoh, Hannes Müller Schmied, Niko Wanders,
# James S. Famiglietti, Johan Rockström, Matti Kummu
#
# Published in Nature Water
# https://doi.org/10.1038/s44221-024-00208-7
# 
# Data availability:
# https://doi.org/10.5281/zenodo.10531807
# 
# Code availability:
# https://github.com/vvirkki/freshwater-pb
# 
# Corresponding authors of the article:
# Miina Porkka (miina.porkka@uef.fi)
# Vili Virkki (vili.virkki@aalto.fi)
# 
# Script author:
# Vili Virkki (vili.virkki@aalto.fi)
#
# Date:
# 15.1.2024
# 
# ------------------------------------------------------------------------------

options(dplyr.summarise.inform = FALSE)

# fill value is mean of all non-NA values
.fill_NA_values <- function(x) {
  xvals <- x[3:length(x)]
  xvals <- ifelse(is.na(xvals), mean(xvals, na.rm = TRUE), xvals)
  ret <- c(x[1], x[2], xvals) %>%
    setNames(names(x))
  return (ret)
}

.detect_cpt <- function(x) {
  xvals <- x[3:length(x)]
  if (var(xvals) < 1e-6) {
    cp <- length(xvals) # minor variance; changepoint set to the end of data
  } else {
    cp <- cpt.meanvar(xvals, method = "AMOC", class = FALSE)["cpt"]
  }
  ret <- c(x[1], x[2], cp) %>%
    setNames(c("x", "y", "cpt"))
  return (ret)
}

# Daily ISIMIP NetCDF -> monthly mean rast -------------------------------------
.daily_to_monthly <- function(file_in, nm_prefix) {

  r <- rast(file_in)
  names(r) <- time(r, format = "days") %>%
    str_replace_all("-", ".") # no dashes to column names
  nm <- names(r)
  yrs <- nm %>%
    substr(1, 4) %>%
    as.numeric() %>%
    unique()

  mnths <- str_pad(1:12, 2, "left", 0)

  collect <- c()
  for (y in 1:length(yrs)) {
    for (m in 1:length(mnths)) {
      sel <- paste0(yrs[y], ".", mnths[m])
      sel_layers <- nm[grepl(sel, nm, fixed = TRUE) & !grepl(".02.29", nm)] # exclude leap days
      month_values <- r %>%
        subset(sel_layers) %>%
        mean() %>%
        setNames(paste0(nm_prefix, sel))
      collect <- c(collect, month_values)
    }
  }
  ret <- rast(collect)
  return (ret)

}

# 10-year daily rasters -> monthly full-timeseries cell data frames ------------
.prepare_df <- function(cfg, period, impactmodel, forcing) {

  # setup
  period_no <- which(cfg$period_label == period)
  ystart <- cfg$period_begin[period_no]
  yend <- cfg$period_end[period_no]

  raw_files <- paste0(cfg$raw_data_root, "/", cfg$variable, "/", period, "/raw") %>%
    list.files(full.names = TRUE)
  files_in <- raw_files[grepl(impactmodel, raw_files)]
  files_in <- files_in[grepl(forcing, files_in)]

  if (all(grepl("daily", files_in))) { # daily data

    data_chunks <- files_in %>%
      lapply(., .daily_to_monthly, nm_prefix = cfg$variable)
    data_rast <- rast(data_chunks)

  } else { # monthly data

    data_rast <- rast(files_in)
    # ISIMIP rootmoist data does not have an adequate time attribute to set names
    nm <- paste0(cfg$variable,
                 sort(rep(seq(ystart, yend, 1), 12)), ".",
                 str_pad(seq(1,12,1), 2, side = "left", pad = "0"))
    names(data_rast) <- nm

  }

  cell_areas <- data_rast[[1]] %>%
    cellSize(unit = "km") %>%
    setNames(c("cellArea_km2")) %>%
    as.data.frame(xy = TRUE)

  impactmodel_out <- str_replace_all(impactmodel, "-", "_")
  forcing_out <- str_replace_all(forcing, "-", "_")
  mids <- str_pad(seq(1,12,1), 2, side = "left", pad = "0")
  data_nm <- names(data_rast)

  for (m in 1:length(mids)) {

    sel <- paste0(".", mids[m])
    sel_layers <- data_nm[grepl(sel, data_nm, fixed = TRUE)]

    # rast to data frame
    mnth_df <- data_rast %>%
      subset(sel_layers) %>%
      as.data.frame(xy = TRUE, cells = TRUE) %>%
      mutate(impactmodel = impactmodel_out,
             forcing = forcing_out,
             intercalibrated = FALSE) %>%
      left_join(cell_areas, by = c("x" = "x", "y" = "y"))

    # check for NA values, fill with mean of non-NA values and write to log if any
    na_values <- mnth_df %>%
      filter(if_any(starts_with(cfg$variable), ~ is.na(.x)))
    if (nrow(na_values) > 0) {

      nvals_mat <- na_values %>%
        select(starts_with(cfg$variable)) %>%
        as.matrix()

      na_filled_values <- mnth_df %>%
        filter(cell %in% na_values$cell) %>%
        select(x, y, starts_with(cfg$variable)) %>%
        as.matrix() %>%
        apply(MARGIN = 1, FUN = .fill_NA_values, simplify = FALSE) %>%
        bind_rows()

      non_filled_values <- mnth_df %>%
        filter(!cell %in% na_values$cell) %>%
        select(x, y, starts_with(cfg$variable))

      mnth_df <- mnth_df %>%
        select(-starts_with(cfg$variable)) %>%
        left_join(bind_rows(na_filled_values, non_filled_values),
                  by = c("x" = "x", "y" = "y"))

      write(paste0("NA monthly ", cfg$variable, " in ", nrow(na_values), " cells ",
                   "in ", length(which(is.na(nvals_mat))), " values"),
            cfg$log_file, append = TRUE)

    }

    # detect changepoints
    cpts <- mnth_df %>%
      select(x, y, starts_with(cfg$variable)) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .detect_cpt, simplify = FALSE) %>%
      bind_rows()

    mnth_df <- mnth_df %>%
      left_join(cpts, by = c("x" = "x", "y" = "y")) %>%
      select(impactmodel, forcing, cell, cellArea_km2, x, y,
             intercalibrated, cpt, everything())

    # prepare output
    out_fldr <- paste0("Data/", cfg$variable, "/", period, "/", cfg$variable,
                       "_monthly",
                       ifelse(impactmodel != "",
                              paste0("/", impactmodel_out), ""),
                       paste0("/", forcing_out))
    if (!dir.exists(out_fldr)) {dir.create(out_fldr, recursive = TRUE)}

    out_file <- paste0(cfg$variable, "_monthly_", period, "_", month.name[m], "_",
                       ifelse(!is.na(impactmodel_out),
                              paste0(impactmodel_out, "_"), ""),
                       forcing_out, ".rds")
    out <- paste0(out_fldr, "/", out_file)
    saveRDS(mnth_df, out)
    write(paste0("saved ", out), cfg$log_file, append = TRUE)

  }
}

# public function
prepare_monthly_data <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("preparing data...", cfg$log_file, append = TRUE)
  periods <- cfg$period_label
  impactmodels <- cfg$impactmodel
  forcings <- cfg$forcing

  for (i in 1:length(periods)) {
    for (j in 1:length(impactmodels)) {
      for (k in 1:length(forcings)) {
        .prepare_df(cfg, periods[i], impactmodels[j], forcings[k])
      }
    }
  }
}
