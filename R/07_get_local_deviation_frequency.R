# ------------------------------------------------------------------------------
# 
# 07_get_local_deviation frequency
# 
# Compute the frequency of local deviations (months with deviations out of all
# months: Fig. 1e) according to config: years described by
# localfreq_period_begin & localfreq_period_end; months by localfreq_season.
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

.check_dpts <- function(x, bound) {

  xvals <- x[3:length(x)]
  if (bound == "dry") {
    ret <- c(x[1], x[2], sum(xvals < 0, na.rm = TRUE) / length(xvals)) %>%
      setNames(c("x", "y", "prc_mnths_blw"))
  } else if (bound == "wet") {
    ret <- c(x[1], x[2], sum(xvals > 0, na.rm = TRUE) / length(xvals)) %>%
      setNames(c("x", "y", "prc_mnths_abv"))
  } else {
    stop("'bound' for '.check_dpts()' should be 'dry' or 'wet'")
  }
  return (ret)

}

.compute_month_shares <- function(cfg, file_departures, ice_df) {

  i_out <- str_replace_all(cfg$impactmodel, "-", "_")
  f_out <- str_replace_all(cfg$forcing, "-", "_")

  dpt <- readRDS(file_departures) %>%
    left_join(ice_df, by = c("x", "y")) %>%
    mutate(isIce = !is.na(isIce))

  ystart <- cfg$localfreq_period_begin[str_detect(file_departures, cfg$period_label)]
  yend <- cfg$localfreq_period_end[str_detect(file_departures, cfg$period_label)]

  sel_years <- paste0("dpt_", cfg$variable, seq(ystart, yend, 1)) %>%
    paste(collapse = "|")

  dryDpt_frequency <- dpt %>%
    select(x, y, matches(sel_years)) %>%
    as.matrix() %>%
    apply(MARGIN = 1, FUN = .check_dpts, bound = "dry", simplify = FALSE) %>%
    bind_rows()

  wetDpt_frequency <- dpt %>%
    select(x, y, matches(sel_years)) %>%
    as.matrix() %>%
    apply(MARGIN = 1, FUN = .check_dpts, bound = "wet", simplify = FALSE) %>%
    bind_rows()

  dpt <- dpt %>%
    select(-starts_with(paste0("dpt_", cfg$variable))) %>%
    left_join(dryDpt_frequency, by = c("x", "y")) %>%
    left_join(wetDpt_frequency, by = c("x", "y"))

  # meta
  cells_dpt <- dpt %>%
    mutate(variable = cfg$variable,
           period = cfg$period_label[str_detect(file_departures, cfg$period_label)],
           ystart = ystart,
           yend = yend) %>%
    select(variable, impactmodel, forcing, period, ystart, yend, x, y,
           prc_mnths_blw, prc_mnths_abv, isIce)

  return (cells_dpt)

}

# public function
get_local_deviation_frequency <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("computing the frequency of local deviations...",
        cfg$log_file, append = TRUE)
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_")

  # check that enough directories exist
  for (p in 1:length(cfg$period_label)) {
    dpt_fldr <- paste0("Data/", cfg$variable, "/", cfg$period_label[p],
                       "/", cfg$variable, "_local_departure_frequency/")
    for (i in 1:length(i_out)) {
      for (j in 1:length(f_out)) {
        outdir <- paste0(dpt_fldr,
                         ifelse(i_out[i] == "", "", paste0(i_out[i], "/")),
                         f_out[j])
        if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
      }
    }
  }

  # query files to process
  i_out <- i_out %>% paste(collapse = "|")
  f_out <- f_out %>% paste(collapse = "|")
  dpt_files <- list.files(paste0("Data/", cfg$variable),
                          recursive = TRUE, full.names = TRUE)
  dpt_files <- dpt_files[grepl(paste0(cfg$variable, "_departures"), dpt_files) &
                         grepl(f_out, dpt_files) &
                         grepl(i_out, dpt_files)]

  # compute and save cell-wise local deviation frequencies
  for (p in 1:length(cfg$period_label)) {

    # ice areas to be excluded
    ice_path <- paste0("Data/", cfg$variable, "/", cfg$period_label[p], "/",
                       cfg$variable, "_departures/", cfg$variable,
                       "_ice_cells_", cfg$period_label[p], "_", cfg$grid_label,
                       ".rds")
    if (!file.exists(ice_path)) {stop("Ice has not been built for this config")}
    ice_df <- readRDS(ice_path)

    for (m in 1:12) {

      month_files <- dpt_files[grepl(month.name[m], dpt_files) &
                               grepl(cfg$period_label[p], dpt_files)]
      month_localfreq <- month_files %>%
        lapply(FUN = .compute_month_shares, cfg = cfg, ice_df = ice_df)

      out_files <- month_files %>%
        lapply(FUN = function(x, cfg) {
          str_replace_all(x, paste0(cfg$variable, "_departures"),
                          paste0(cfg$variable, "_local_departure_frequency"))},
          cfg = cfg) %>%
        lapply(FUN = function(x, cfg) {
          str_replace_all(x, paste0("local_departure_frequency_", cfg$period_label[p]),
                          paste0("local_departure_frequency_", cfg$period_label[p],
                                 "_", cfg$localfreq_period_begin[p],
                                 "_", cfg$localfreq_period_end[p]))},
          cfg = cfg)

      for (i in 1:length(month_localfreq)) {
        saveRDS(month_localfreq[[i]], out_files[[i]])
        write(paste0("saved ", out_files[[i]]), cfg$log_file, append = TRUE)
      }

    }
  }
}
