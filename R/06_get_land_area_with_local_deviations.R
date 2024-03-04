# ------------------------------------------------------------------------------
# 
# 06_get_land_area_with_local_deviations
# 
# Aggregate local deviations globally or within areas prescribed in config and
# given as tif raster in Data/ardiv to get percentage of land area with local
# deviations (Fig. 1c).
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

.classify_dpt <- function(dpt) { ifelse(dpt < 0, "below", "above") }

.aggregate_dpts <- function(cfg, dpts_df, ice_df, ardiv_id, dpts_period,
                            ardiv_df = NULL) {

  # exclusion of ice cells
  dpts_df <- dpts_df %>%
    left_join(ice_df, by = c("x", "y")) %>%
    filter(is.na(isIce)) %>%
    select(-isIce)

  # filtering of cells within a specific area
  if (!is.null(ardiv_df)) {
    dpts_df <- dpts_df %>%
      left_join(ardiv_df, by = c("x", "y")) %>%
      filter(id == ardiv_id)
    area_meta <- paste0(cfg$areal_division, "_id", ardiv_id)
  } else {
    area_meta <- "global"
  }

  # binary departures
  dpt_classified <- dpts_df %>%
    mutate(across(starts_with("dpt_"), ~ .classify_dpt(.x)))

  # land area shares
  area_dpt <- dpt_classified %>%
    select(x, y, cellArea_km2, starts_with("dpt_")) %>%
    pivot_longer(-c(x, y, cellArea_km2),
                 names_to = "timestep", values_to = "dptClass") %>%
    group_by(timestep, dptClass) %>%
    summarise(areaInClass = sum(cellArea_km2)) %>%
    mutate(landShareInClass = areaInClass / sum(areaInClass),
           dptClass = factor(dptClass, levels = c("above", "below", NA))) %>%
    ungroup()

  # ensure all dpt classes being present in data
  area_dpt_full <- area_dpt %>%
    expand(timestep, dptClass) %>%
    left_join(area_dpt, by = c("timestep" = "timestep",
                               "dptClass" = "dptClass")) %>%
    mutate(landShareInClass = ifelse(is.na(landShareInClass), 0, landShareInClass))

  # meta
  area_dpt_full <- area_dpt_full %>%
    mutate(variable = cfg$variable,
           timestep = str_replace(timestep, paste0("dpt_", cfg$variable), ""),
           period = dpts_period,
           impactmodel = unique(dpts_df$impactmodel),
           forcing = unique(dpts_df$forcing),
           area = area_meta) %>%
    tidyr::separate(timestep, into = c("year", "month")) %>%
    mutate(year = as.numeric(year),
           month = as.numeric(month)) %>%
    arrange(year) %>%
    select(-areaInClass) %>%
    select(variable, impactmodel, forcing, period, area, everything())

  return (area_dpt_full)

}

# public function
get_land_area_with_local_deviations <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("computing percentages of land area with local deviations...",
        cfg$log_file, append = TRUE)
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_")

  # run for a given areal division or global if nothing is given
  ardiv_lbl <- ifelse(cfg$areal_division == "", "global", cfg$areal_division)
  if (ardiv_lbl != "global") {
    # areal division must be given as 0.5-degree tif raster with one value = one area
    ardiv_df <- rast(paste0("data/ardiv/", cfg$areal_division, ".tif")) %>%
      as.data.frame(xy = TRUE) %>%
      as_tibble()
    ardiv_ids <- unique(ardiv_df$id)
  } else {
    ardiv_df <- NULL
    ardiv_ids <- c(NA)
  }
  write(paste0("processing areal division ", ardiv_lbl, "..."),
        cfg$log_file, append = TRUE)

  # check that enough directories exist
  for (p in 1:length(cfg$period_label)) {
    dpt_fldr <- paste0("Data/", cfg$variable, "/", cfg$period_label[p],
                       "/", cfg$variable, "_", ardiv_lbl, "_departure_aggregates/")
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

  # compute and save percentages of land area with local deviations
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
      month_data <- month_files %>%
        lapply(FUN = readRDS)

      ardiv_aggregates <- vector(mode = "list", length = length(ardiv_ids))
      for (a in 1:length(ardiv_ids)) {
        month_aggregates <- month_data %>%
          lapply(FUN = .aggregate_dpts, cfg = cfg, ice_df = ice_df,
                 ardiv_id = ardiv_ids[a], dpts_period = cfg$period_label[p],
                 ardiv_df = ardiv_df)
        ardiv_aggregates[[a]] <- month_aggregates
      }

      out_id <- ifelse(!cfg$areal_division == "", cfg$areal_division, "global")
      out_files <- month_files %>%
        lapply(FUN = function(x, cfg) {
          str_replace_all(x, paste0(cfg$variable, "_departures"),
                          paste0(cfg$variable, "_", out_id, "_departure_aggregates"))},
          cfg = cfg)

      for (i in 1:length(month_files)) {

        all_ardiv_aggregates <- vector(mode = "list", length = length(ardiv_aggregates))
        for (j in 1:length(ardiv_aggregates)) {
          all_ardiv_aggregates[[j]] <- ardiv_aggregates[[j]][[i]]
        }

        # dimensions should be 3 classes x n years per period
        data_coverage <- 3 * (cfg$period_end[p] - cfg$period_begin[p] + 1)
        check_coverage <- all_ardiv_aggregates %>%
          lapply(nrow) %>%
          unlist()
        not_covered <- which(check_coverage != data_coverage)

        # log areas which are not covered by departures data
        if (length(not_covered) > 0) {
          write(paste0("data not covering area id ", ardiv_ids[not_covered]),
                cfg$log_file, append = TRUE)
        }
        all_ardiv_aggregates <- all_ardiv_aggregates[check_coverage == data_coverage] %>%
          bind_rows()
        saveRDS(all_ardiv_aggregates, out_files[[i]])
        write(paste0("saved ", out_files[[i]]), cfg$log_file, append = TRUE)

      }
    }
  }
}
