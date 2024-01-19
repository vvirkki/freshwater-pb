# ------------------------------------------------------------------------------
# 
# output02_local_deviation_frequency_change
# 
# Compare the frequency of local deviations between baseline and comparison
# periods described in config, test significance of differences, output grid
# cell classifications and class area summaries (plot, csv).
#
# Figures using these outputs: 3, 5
# Extended Data Figures using these outputs: 5, 6, 7
# 
# Analysis script for publication:
# 
# "Global water cycle shifts far beyond pre-industrial conditions - planetary
# boundary for freshwater change transgressed"
#
# Miina Porkka, Vili Virkki, Lan Wang-Erlandsson, Dieter Gerten, Tom Gleeson,
# Chinchu Mohan, Ingo Fetzer, Fernando Jaramillo, Arie Staal, Sofie te Wierik,
# Arne Tobian, Ruud van der Ent, Petra Döll, Martina Flörke, Simon N. Gosling,
# Naota Hanasaki, Yusuke Satoh, Hannes Müller Schmied, Niko Wanders,
# James S. Famiglietti, Johan Rockström, Matti Kummu
#
# Published in Nature Water
# DOI: TO-ADD
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

.iceId <- -999
.exclusions <- c(.iceId)
.map_varVoc <- c("streamflow", "root-zone soil moisture") %>%
  setNames(c("dis", "rootmoist"))

.test_proportions <- function(x, ntot, pval, bound) {

  if (bound %in% c("blw", "abv")) {
    b <- ifelse(bound == "blw", "dryExit", "wetExit")
    p <- ifelse(bound == "blw", "a", "b")
  } else {
    stop ("'bound' should be 'blw' or 'abv'")
  }

  if (any(x[3:4] < 0) | any(is.na(x[3:4]))) {
    changeclass <- NA_character_
  } else {
    suppressWarnings({
      # H1: proportion during comparison period is significantly greater than during baseline period
      h0_baseline_prop_greater <- prop.test(x[3:4] * ntot, ntot, alternative = "less")
      # H1: proportion during comparison period is significantly less than during baseline period
      h0_baseline_prop_less <- prop.test(x[3:4] * ntot, ntot, alternative = "greater")
    })
    changeclass <- case_when(
      # H1 of this test stands if p.value < pval --> significant increase in proportion
      h0_baseline_prop_greater$p.value < pval ~ paste0(p, "a_", b, "Increase"),
      # H1 of this test stands if p.value < pval --> significant decrease in proportion
      h0_baseline_prop_less$p.value < pval ~ paste0(p, "b_", b, "Decrease"),
      # both tests pass --> no significant change in proportion
      TRUE ~ paste0(p, "c_", "no", b, "Change")
    )
  }

  ret <- c(x[1], x[2], changeclass) %>% setNames(c("x", "y", "changeclass"))
  return (ret)

}

.prepare_for_map <- function(cfg, mapdata, test_changes = FALSE) {

  # median over ensemble members
  mapdata_medians <- mapdata %>%
    group_by(period, x, y, isIce) %>%
    summarise(median_prc_blw = median(prc_mnths_blw),
              median_prc_abv = median(prc_mnths_abv)) %>%
    ungroup() %>%
    mutate(across(matches("prc_blw|prc_abv"),
           ~ case_when(
             isIce ~ .iceId,
             TRUE ~ .x
           )))

  if (test_changes) {

    nmtot <- mapdata %>%
      select(period, nmtot_sum) %>%
      distinct()
    ntot <- c(nmtot %>% filter(period == cfg$baseline_period) %>% pull(nmtot_sum),
              nmtot %>% filter(!period == cfg$baseline_period) %>% pull(nmtot_sum))

    stat_changes_prep <- mapdata_medians %>%
      mutate(period = ifelse(period == cfg$baseline_period,
                             "baseline", "comparison")) %>%
      pivot_wider(id_cols = c(x, y, isIce),
                  names_from = period,
                  values_from = matches("prc"))

    write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
    write("testing dry local deviation frequency change significance...",
          cfg$log_file, append = TRUE)

    blwChange <- stat_changes_prep %>%
      select(x, y, median_prc_blw_baseline, median_prc_blw_comparison) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .test_proportions, ntot = ntot,
            pval = cfg$localfreq_chg_pval, bound = "blw", simplify = FALSE) %>%
      bind_rows() %>%
      rename(blwChange = changeclass) %>%
      mutate(x = as.numeric(x),
             y = as.numeric(y))

    write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
    write("testing wet local deviation frequency change significance...",
          cfg$log_file, append = TRUE)

    abvChange <- stat_changes_prep %>%
      select(x, y, median_prc_abv_baseline, median_prc_abv_comparison) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .test_proportions, ntot = ntot,
            pval = cfg$localfreq_chg_pval, bound = "abv", simplify = FALSE) %>%
      bind_rows() %>%
      rename(abvChange = changeclass) %>%
      mutate(x = as.numeric(x),
             y = as.numeric(y))

    stat_changes <- stat_changes_prep %>%
      select(x, y, isIce) %>%
      left_join(blwChange, by = c("x" = "x", "y" = "y")) %>%
      left_join(abvChange, by = c("x" = "x", "y" = "y"))

  }

  baseline_shares <- mapdata_medians %>%
    filter(period == cfg$period_label[which(cfg$period_label == cfg$baseline_period)]) %>%
    select(-period) %>%
    rename_with(~ paste0("baseline_", .), -c(x, y))

  comparison_shares <- mapdata_medians %>%
    filter(period ==  cfg$period_label[which(cfg$period_label != cfg$baseline_period)]) %>%
    select(-period) %>%
    rename_with(~ paste0("comparison_", .), -c(x,y))

  difference_shares <- baseline_shares %>%
    left_join(comparison_shares, by = c("x", "y")) %>%
    mutate(diff_isIce = baseline_isIce | comparison_isIce,
           diff_median_prc_blw = comparison_median_prc_blw - baseline_median_prc_blw,
           diff_median_prc_abv = comparison_median_prc_abv - baseline_median_prc_abv) %>%
    select(x, y, starts_with("diff")) %>%
    mutate(across(matches("prc_blw|prc_abv"),
                  ~ case_when(
                    diff_isIce ~ .iceId,
                    TRUE ~ .x
                  )))

  ret <- list(baseline_shares, comparison_shares, difference_shares) %>%
    lapply(function(x) { x %>% select(-matches("isIce")) }) %>%
    setNames(c("baseline", "comparison", "difference"))

  if (test_changes) {
    ret$stat_changes <- stat_changes
  } else {
    ret$stat_changes <- NA
  }
  return (ret)

}

.tbl_to_rast <- function(tbl_in, variable) {

  r <- tbl_in %>%
    select(x, y, !!variable) %>%
    as.matrix() %>%
    rast(type = "xyz", crs = "+proj=longlat")
  return (r)

}

.do_changeclass_map <- function(cfg, visopts, chg_data, season_abbr, chgmagn) {

  if (chgmagn) {
    b <- "#4E0F2B"
    map_ice_col <- "grey80"
    if (cfg$variable == "dis") {
      map_colours <- c(b, b, "#B7AC07", b, b, "#FAF28D", "#166673", "#a1ced6", "grey90")
    } else {
      map_colours <- c(b, b, "#806227", b, b, "#c6b496", "#0F6D37", "#9fbdae", "grey90")
    }
    blw_map <- c("aa_dryExitIncrease_a_major", "aa_dryExitIncrease_b_minor",
                 "ca_noFrequencyChange_or_FrequencyDecrease")
    abv_map <- c("ba_wetExitIncrease_a_major", "ba_wetExitIncrease_b_minor",
                 "ca_noFrequencyChange_or_FrequencyDecrease")
  } else {
    map_colours <- c("#000000", "#f3b300", "#b36600", "#509dc2", "#f3f3f3",
    "#dfeeef", "#376387", "#f4ecdc", "#cccccc")
    map_ice_col <- "grey60"
    blw_map <- c("aa_dryExitIncrease", "ab_dryExitDecrease", "ac_nodryExitChange")
    abv_map <- c("ba_wetExitIncrease", "bb_wetExitDecrease", "bc_nowetExitChange")
  }

  possible_class_combinations <- tidyr::expand_grid(blw_map, abv_map) %>%
    mutate(class_levels = paste0(blw_map, "_", abv_map))
  possible_class_map <- possible_class_combinations %>%
    pull(class_levels)

  # map categories to numeric classes
  class_map <- 1:length(unique(chg_data$class)) %>%
    setNames(chg_data$class %>% as.factor() %>% levels())

  vis_data <- chg_data %>%
    mutate(class_numeric = case_when(
      isIce ~ .iceId,
      TRUE ~ as.numeric(class_map[class]))) %>%
    filter(class != "NA_NA" | class_numeric %in% .exclusions)

  # do raster
  rast_classes <- .tbl_to_rast(vis_data, "class_numeric")
  rast_area <- cellSize(rast_classes, unit = "km") %>%
    as.data.frame(xy = TRUE) %>%
    rename(area_km2 = 3)

  # compute areas in each class
  class_areas <- vis_data %>%
    left_join(rast_area, by = c("x", "y")) %>%
    filter(class_numeric != .iceId) %>% # ice-free land area only
    mutate(class = factor(class, levels = possible_class_map)) %>%
    group_by(class, .drop = FALSE) %>%
    summarise(area_tot_km2 = sum(area_km2)) %>%
    mutate(area_share_perc = 100 * (area_tot_km2 / sum(.$area_tot_km2)))

  # labels
  lbl_included <- class_areas %>%
    filter(class %in% unique(vis_data$class)) %>%
    mutate(class = as.character(class))

  lbl_shr_areas <- lbl_included$class %>%
    str_replace("aa_|ab_|ac_", "") %>%
    str_replace("_ba_|_bb_|_bc_", " ") %>%
    str_replace("_ca_", " ") %>%
    str_replace("ca_", "") %>%
    str_replace_all("Exit", "DeviationFrequency") %>%
    str_replace_all("_a_|_b_", " ") %>%
    paste0(., " (", round(lbl_included$area_share_perc, 3), "%)")

  lbl <- c("ice")
  lbl <- c(lbl[.exclusions %in% values(rast_classes)], lbl_shr_areas)

  # colors
  pal <- c(map_ice_col)
  pal <- pal[.exclusions %in% values(rast_classes)]
  names(map_colours) <- possible_class_map
  pal <- c(pal, map_colours[names(map_colours) %in% unique(vis_data$class)] %>% as.character())

  # map title
  ib <- which(cfg$period_label == cfg$baseline_period)
  ic <- which(cfg$period_label != cfg$baseline_period)
  map_title <- paste0("Change in the frequency of local ", .map_varVoc[cfg$variable],
                      " deviations (", cfg$localfreq_period_begin[ic], "-",
                      cfg$localfreq_period_end[ic], " compared to ",
                      cfg$localfreq_period_begin[ib], "-", cfg$localfreq_period_end[ib],
                      "; months ", season_abbr, ")")

  rast_classes_rob <- rast_classes %>%
    project("+proj=robin", mask = TRUE, method = "near")

  mm <- tm_shape(rast_classes_rob) +
    tm_raster(title = "Changeclass category (% of ice-free land area)",
              style = "cat",
              labels = lbl,
              palette = pal) +
    tm_layout(main.title = map_title,
              main.title.size = visopts$map_title_size,
              frame = FALSE,
              legend.outside = TRUE,
              legend.outside.position = "bottom",
              legend.outside.size = visopts$map_legend_size)

  return (list(mm, class_areas, possible_class_combinations) %>%
            setNames(c("map", "class_areas", "class_map")))

}

# public function
output_local_deviation_frequency_change <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("outputting local deviation frequency changes...",
        cfg$log_file, append = TRUE)

  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")

  visopts_cfg <- paste0("configs/", cfg$visopts_cfg)
  if (!file.exists(visopts_cfg)) {
    stop("visopts not defined for this config")
  } else {
    visopts <- .read_visopts(visopts_cfg)
  }

  # get months to be included in maps
  if (!any(cfg$localfreq_season == "")) {
    season_months <- cfg$localfreq_season # months specified by config
  } else {
    season_months <- month.name # all months
  }

  season_abbr <- season_months %>%
    lapply(function(x) {str_sub(x, 1, 1)}) %>%
    unlist() %>%
    paste(collapse = "")
  localfreq_period_abbr <- paste0(cfg$localfreq_period_begin, "_", cfg$localfreq_period_end)

  tstamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_fldr <- paste0("output/", cfg$variable, "_local_deviation_frequency_changes_",
                        paste(localfreq_period_abbr, collapse = "_to_"), "_", season_abbr, "_",
                        paste(cfg$impactmodel, collapse = "_"), "_", tstamp)
  if (!dir.exists(output_fldr)) { dir.create(output_fldr, recursive = TRUE) }

  # query files to process
  mapdata_files <- list.files(paste0("Data/", cfg$variable),
                              recursive = TRUE, full.names = TRUE)
  map_period <- localfreq_period_abbr %>%
    paste(collapse = "|")
  mapdata_files <- mapdata_files[grepl("local_departure_frequency", mapdata_files) &
                                 grepl(map_period, mapdata_files) &
                                 grepl(f_out, mapdata_files) &
                                 grepl(i_out, mapdata_files)]

  # collect monthly data into one
  mapdata_season_months <- tibble()
  for (m in 1:12) {

    mapdata_month <- mapdata_files[grepl(month.name[m], mapdata_files)] %>%
      lapply(readRDS) %>%
      bind_rows() %>%
      mutate(month = month.name[m])

    mapdata_season_months <- mapdata_season_months %>%
      bind_rows(mapdata_month)

  }
  nensmem <- mapdata_season_months %>%
    group_by(impactmodel, forcing) %>%
    n_groups()

  # share of departing months out of all months for each cell and ensemble member
  mapdata_season_months_shares <- mapdata_season_months %>%
    filter(month %in% season_months) %>%
    mutate(nmblw = prc_mnths_blw * (yend - ystart + 1),
           nmabv = prc_mnths_abv * (yend - ystart + 1),
           nmtot = (yend - ystart + 1)) %>%
    group_by(x, y, impactmodel, forcing, period, isIce) %>%
    summarise(nmblw_sum = sum(nmblw),
              nmabv_sum = sum(nmabv),
              nmtot_sum = sum(nmtot)) %>%
    ungroup() %>%
    mutate(prc_mnths_blw = nmblw_sum / nmtot_sum,
           prc_mnths_abv = nmabv_sum / nmtot_sum)

  # prepare and save departure frequencies
  vis_mapdata_season_months <- mapdata_season_months_shares %>%
    .prepare_for_map(cfg = cfg, test_changes = TRUE)

  # statistically significant changes (all changes; no magnitude)
  chgmap_data <- vis_mapdata_season_months$stat_changes %>%
    mutate(class = paste0(blwChange, "_", abvChange))

  chgmap <- .do_changeclass_map(cfg, visopts, chgmap_data, season_abbr,
                                chgmagn = FALSE)

  # map pdf output
  chgmap_out <- paste0(output_fldr, "/", cfg$variable,
                       "_changeclasses_all_frequency_changes_ensn",
                       nensmem, ".pdf")
  tmap_save(chgmap$map, chgmap_out,
            width = 1.5 * visopts$map_out_width,
            height = 1.5 * visopts$map_out_height,
            units = visopts$map_out_units)
  write(paste0("saved ", chgmap_out, "..."), cfg$log_file, append = TRUE)

  # csv output
  chgmap_data_out <- chgmap_data %>%
    mutate(blwChange = str_replace_all(blwChange, "aa_|ab_|ac_|dry", ""),
           blwChange = str_replace_all(blwChange, "Exit", "Frequency"),
           abvChange = str_replace_all(abvChange, "ba_|bb_|bc_|wet", ""),
           abvChange = str_replace_all(abvChange, "Exit", "Frequency")) %>%
    rename(dry_local_deviations = blwChange,
           wet_local_deviations = abvChange) %>%
    select(-class)

  chgmap_data_out_path <- paste0(output_fldr, "/", cfg$variable,
                                 "_changeclasses_in_cells_all_frequency_changes_ensn",
                                 nensmem, ".csv")
  write.csv(chgmap_data_out, chgmap_data_out_path, row.names = FALSE)
  write(paste0("saved ", chgmap_data_out_path, "..."), cfg$log_file, append = TRUE)

  # csv output
  chgclass_areas <- chgmap$class_areas %>%
    left_join(chgmap_data %>% select(blwChange, abvChange, class) %>% distinct(),
              by = "class") %>%
    mutate(blwChange = str_replace_all(blwChange, "aa_|ab_|ac_|dry", ""),
           blwChange = str_replace_all(blwChange, "Exit", "Frequency"),
           abvChange = str_replace_all(abvChange, "ba_|bb_|bc_|wet", ""),
           abvChange = str_replace_all(abvChange, "Exit", "Frequency")) %>%
    rename(dry_local_deviations = blwChange,
           wet_local_deviations = abvChange) %>%
    select(dry_local_deviations, wet_local_deviations, area_tot_km2, area_share_perc)

  chgclass_areas_out <- paste0(output_fldr, "/", cfg$variable,
                               "_changeclass_areas_all_frequency_changes_ensn",
                               nensmem, ".csv")
  write.csv(chgclass_areas, chgclass_areas_out, row.names = FALSE)
  write(paste0("saved ", chgclass_areas_out, "..."), cfg$log_file, append = TRUE)

  # statistically significant changes (increases only; with magnitude)
  change_magnitude <- vis_mapdata_season_months$difference %>%
    mutate(across(-c(x,y), ~ ifelse(abs(.x) > 0.05, "a_major", "b_minor"))) %>%
    setNames(c("x", "y", "blwMagn", "abvMagn"))

  chgmap_data_magnitudes <- chgmap_data %>%
    left_join(change_magnitude, by = c("x", "y")) %>%
    mutate(across(ends_with("Change"),
                  ~ ifelse(grepl("Increase", .x), .x, "ca_noFrequencyChange_or_FrequencyDecrease")),
           blwChange = ifelse(grepl("Increase", blwChange),
                              paste0(blwChange, "_", blwMagn),
                              blwChange),
           abvChange = ifelse(grepl("Increase", abvChange),
                              paste0(abvChange, "_", abvMagn),
                              abvChange)) %>%
    mutate(class = paste0(blwChange, "_", abvChange)) %>%
    select(x, y, isIce, blwChange, abvChange, blwMagn, abvMagn, class)

  chgmap_magn <- .do_changeclass_map(cfg, visopts, chgmap_data_magnitudes, season_abbr,
                                     chgmagn = TRUE)

  # map pdf output
  chgmap_magn_out <- paste0(output_fldr, "/", cfg$variable,
                            "_changeclasses_frequency_increases_only_with_magnitude_ensn",
                            nensmem, ".pdf")
  tmap_save(chgmap_magn$map, chgmap_magn_out,
            width = 1.5 * visopts$map_out_width,
            height = 1.5 * visopts$map_out_height,
            units = visopts$map_out_units)
  write(paste0("saved ", chgmap_magn_out, "..."), cfg$log_file, append = TRUE)

  # csv output
  chgmap_data_magn_out <- chgmap_data_magnitudes %>%
    mutate(blwChange = str_replace_all(blwChange, "aa_|ab_|ac_|ca_|dry|_a_major|_b_minor", ""),
           blwChange = str_replace_all(blwChange, "Exit", "Frequency"),
           abvChange = str_replace_all(abvChange, "ba_|bb_|bc_|ca_|wet|_a_major|_b_minor", ""),
           abvChange = str_replace_all(abvChange, "Exit", "Frequency"),
           blwMagn = str_replace_all(blwMagn, "a_|b_", ""),
           blwMagn = ifelse(grepl("noFrequencyChange", blwChange), NA, blwMagn),
           abvMagn = str_replace_all(abvMagn, "a_|b_", ""),
           abvMagn = ifelse(grepl("noFrequencyChange", abvChange), NA, abvMagn),
           across(-c(x, y, isIce), ~ ifelse(isIce, NA, .x))) %>%
    rename(dry_local_deviations = blwChange,
           dry_frequency_increase_magnitude = blwMagn,
           wet_local_deviations = abvChange,
           wet_frequency_increase_magnitude = abvMagn) %>%
    select(x, y, isIce, starts_with("dry"), starts_with("wet"))

  chgmap_data_magn_out_path <- paste0(output_fldr, "/", cfg$variable,
                                      "_changeclasses_in_cells_frequency_increases_only_with_magnitude_ensn",
                                      nensmem, ".csv")
  write.csv(chgmap_data_magn_out, chgmap_data_magn_out_path, row.names = FALSE)
  write(paste0("saved ", chgmap_data_magn_out_path, "..."), cfg$log_file, append = TRUE)

  # csv output
  chgclass_areas_magn <- chgmap_magn$class_areas %>%
    left_join(chgmap_magn$class_map %>% rename(blwChange = blw_map, abvChange = abv_map),
              by = c("class" = "class_levels")) %>%
    mutate(blwChange = str_replace_all(blwChange, "aa_|ab_|ac_|ca_|dry", ""),
           blwChange = str_replace_all(blwChange, "Exit", "Frequency"),
           abvChange = str_replace_all(abvChange, "ba_|bb_|bc_|ca_|wet", ""),
           abvChange = str_replace_all(abvChange, "Exit", "Frequency")) %>%
    separate_wider_delim(blwChange, delim = regex("_a_|_b_"), names = c("blwChange", "blwMagn"),
                         too_few = "align_start") %>%
    separate_wider_delim(abvChange, delim = regex("_a_|_b_"), names = c("abvChange", "abvMagn"),
                         too_few = "align_start") %>%
    rename(dry_local_deviations = blwChange,
           wet_local_deviations = abvChange,
           dry_frequency_increase_magnitude = blwMagn,
           wet_frequency_increase_magnitude = abvMagn) %>%
    select(starts_with("dry"), starts_with("wet"), area_tot_km2, area_share_perc)

  chgclass_areas_magn_out <- paste0(output_fldr, "/", cfg$variable,
                                    "_changeclass_areas_frequency_increases_only_with_magnitude_ensn",
                                    nensmem, ".csv")
  write.csv(chgclass_areas_magn, chgclass_areas_magn_out, row.names = FALSE)
  write(paste0("saved ", chgclass_areas_magn_out, "..."), cfg$log_file, append = TRUE)

}
