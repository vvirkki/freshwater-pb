# ------------------------------------------------------------------------------
# 
# plot02_combined_local_deviation_frequency_increases
# 
# Combine detected increases in local deviation frequency for streamflow and
# soil moisture and plot classified raster (grid) map.
#
# Figures using these outputs: 5
# Extended Data Figures using these outputs: 7
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

# the latest output folder containing local deviation frequency changes
.find_latest_grid_output <- function(variable, baseline_period, comparison_period, months_abb) {

  opf <- list.files("output", full.names = TRUE)

  latest_fldr <- opf[grepl(variable, opf) &
                     grepl(baseline_period, opf) &
                     grepl(comparison_period, opf) &
                     grepl(months_abb, opf)] %>%
    tibble(fldr = .) %>%
    mutate(ctime = file.info(fldr)$ctime) %>%
    arrange(desc(ctime)) %>%
    slice(1) %>%
    pull(fldr)

  return (latest_fldr)

}


.tbl_to_rast <- function(tbl_in, variable) {

  r <- tbl_in %>%
    select(x, y, !!variable) %>%
    as.matrix() %>%
    rast(type = "xyz", crs = "+proj=longlat")
  return (r)

}

.run_plot02 <- function(bw_localfreq_fldr, gw_localfreq_fldr, yrs,
                         overlay_basins = NULL) {

  comb_map <- c("a_dry_FrequencyIncrease", "b_wet_FrequencyIncrease",
                "a_dry_FrequencyIncrease_b_wet_FrequencyIncrease", "")
  possible_class_combinations <- expand_grid(paste0("a_BW_", comb_map), paste0("b_GW_", comb_map)) %>%
    setNames(c("BW", "GW")) %>%
    mutate(BW = ifelse(BW == "a_BW_", NA, BW),
           GW = ifelse(GW == "b_GW_", NA, GW),
           class = paste0(BW, "_", GW),
           class = str_replace_all(class, "NA_NA|NA_|_NA", ""),
           class = ifelse(class == "", "c_noChanges", class))
  possible_class_map <- possible_class_combinations %>%
    pull(class)

  bwdata <- read_csv(paste0(bw_localfreq_fldr, "/dis_changeclasses_in_cells_frequency_increases_only_with_magnitude_ensn23.csv"),
                     show_col_types = FALSE) %>%
    mutate(bwClass = paste0("a_dry_", dry_local_deviations, "_",
                            "b_wet_", wet_local_deviations),
           bwClass = str_replace_all(bwClass, "noFrequencyChange_or_FrequencyDecrease", "NA"),
           bwClass = str_replace_all(bwClass, "a_dry_NA_|_b_wet_NA|b_wet_NA", ""),
           bwClass = ifelse(bwClass != "", paste0("a_BW_", bwClass), NA)) %>%
    rename(bwIsIce = isIce) %>%
    rename_with(~ paste0("streamflow_", .x), matches("dry_|wet_")) %>%
    select(-matches("magnitude"))

  gwdata <- read_csv(paste0(gw_localfreq_fldr, "/rootmoist_changeclasses_in_cells_frequency_increases_only_with_magnitude_ensn15.csv"),
                     show_col_types = FALSE) %>%
    mutate(gwClass = paste0("a_dry_", dry_local_deviations, "_",
                            "b_wet_", wet_local_deviations),
           gwClass = str_replace_all(gwClass, "noFrequencyChange_or_FrequencyDecrease", "NA"),
           gwClass = str_replace_all(gwClass, "a_dry_NA_|_b_wet_NA|b_wet_NA", ""),
           gwClass = ifelse(gwClass != "", paste0("b_GW_", gwClass), NA)) %>%
    rename(gwIsIce = isIce) %>%
    rename_with(~ paste0("soilmoisture_", .x), matches("dry_|wet_")) %>%
    select(-matches("magnitude"))

  combined_data <- bwdata %>%
    select(x, y, bwIsIce, bwClass) %>%
    left_join(gwdata %>% select(x, y, gwIsIce, gwClass),
              by = c("x", "y")) %>%
    mutate(class = paste0(bwClass, "_", gwClass),
           class = str_replace_all(class, "NA_NA|NA_|_NA", ""),
           class = ifelse(class == "", "c_noChanges", class))

  class_map <- 1:length(unique(combined_data$class)) %>%
    setNames(combined_data$class %>% as.factor() %>% levels())

  vis_data <- combined_data %>%
    mutate(class_numeric = case_when(
      bwIsIce|gwIsIce ~ -999,
      TRUE ~ as.numeric(class_map[class]))
    )

  rast_classes <- .tbl_to_rast(vis_data, "class_numeric")
  rast_area <- cellSize(rast_classes, unit = "km") %>%
    as.data.frame(xy = TRUE) %>%
    rename(area_km2 = 3)

  class_areas <- combined_data %>%
    left_join(rast_area, by = c("x", "y")) %>%
    filter(!bwIsIce & !gwIsIce) %>%
    mutate(class = factor(class, levels = possible_class_map)) %>%
    group_by(class, .drop = FALSE) %>%
    summarise(area_tot_km2 = sum(area_km2)) %>%
    mutate(area_share_perc = 100 * round(area_tot_km2 / sum(.$area_tot_km2), 4)) %>%
    arrange(-area_share_perc) %>%
    left_join(possible_class_combinations, by = "class") %>%
    mutate(BW = str_replace_all(BW, "a_BW_|a_|b_", ""),
           GW = str_replace_all(GW, "b_GW_|a_|b_", ""),
           across(c(BW, GW), ~ str_replace_all(.x, "dry_", "dryDeviation")),
           across(c(BW, GW), ~ str_replace_all(.x, "wet_", "wetDeviation")),
           across(c(BW, GW), ~ str_replace_all(.x, "_", "_and_"))) %>%
    rename(streamflow_local_deviations = BW,
           soilmoisture_local_deviations = GW)

  rast_classes_rob <- rast_classes %>%
    project("+proj=robin", mask = TRUE, method = "near")

  oneChange <- "#c8d6bf" # change only in BW or GW
  otherChng <- "#6b0803" # "other combination"

  colour_map <- c(oneChange, "#b37b1b", "#fcb03d", "#c744a8",
                  oneChange, otherChng, otherChng, "#383bc9",
                  oneChange, otherChng, otherChng, "#00d0db",
                  oneChange, oneChange, oneChange, "#e5e5e5") %>%
    setNames(sort(possible_class_map))

  colour_map <- colour_map[names(colour_map) %in% unique(vis_data$class)]

  lbl <- names(colour_map) %>%
    str_replace_all("BW_", "streamflow_") %>%
    str_replace_all("GW_", "soilmoisture_") %>%
    str_replace_all("dry_", "dryDeviation") %>%
    str_replace_all("wet_", "wetDeviation") %>%
    str_replace_all("a_|b_|c_", "")

  mm <- tm_shape(rast_classes_rob) +
    tm_raster(title = "class",
              style = "cat",
              breaks = c(-999, class_map),
              labels = c("ice", lbl),
              palette = c("grey80", as.character(colour_map))) +
    tm_layout(frame = FALSE,
              legend.outside = TRUE,
              legend.outside.position = "bottom")

  if (!is.null(overlay_basins)) {
    mm <- mm +
      tm_shape(overlay_basins) +
      tm_borders(col = "grey20")
  }

  out_fldr <- "output/combined_local_deviation_frequency_increases/"
  tmap_save(mm, paste0(out_fldr, "combined_frequency_increase_changeclasses_",
                       yrs, ".pdf"),
            width = 1.5 * 85,
            height = 1.5 * 55,
            units = "mm")

  combined_classes_out <- vis_data %>%
    left_join(class_areas %>% select(class,
                                     streamflow_local_deviations,
                                     soilmoisture_local_deviations), by = "class") %>%
    mutate(across(ends_with("local_deviations"),
                  ~ ifelse(is.na(.x), "no_DeviationFrequencyIncrease", .x)),
           isIce = bwIsIce | gwIsIce,
           across(ends_with("local_deviations"),
                  ~ ifelse(isIce, NA, .x))) %>%
    select(x, y, isIce, streamflow_local_deviations, soilmoisture_local_deviations)

  write.csv(combined_classes_out,
            paste0(out_fldr, "combined_frequency_increase_changeclasses_in_cells_",
                   yrs, ".csv"),
            row.names = FALSE)

  write.csv(class_areas %>%
              select(streamflow_local_deviations, soilmoisture_local_deviations,
                     area_tot_km2, area_share_perc) %>%
              arrange(streamflow_local_deviations, soilmoisture_local_deviations) %>%
              mutate(across(ends_with("local_deviations"),
                            ~ ifelse(is.na(.x), "no_DeviationFrequencyIncrease", .x))),
            paste0(out_fldr, "combined_frequency_increase_changeclass_areas_",
                   yrs, ".csv"),
            row.names = FALSE)

}

# public function
plot_combined_local_deviation_frequency_increases <- function() {

  if (!dir.exists("output/combined_local_deviation_frequency_increases")) {
    dir.create("output/combined_local_deviation_frequency_increases", recursive = TRUE)
  }
  all_months <- month.name %>%
    lapply(str_sub, start = 1, end = 1) %>%
    unlist() %>%
    paste(collapse = "")

  bw_1976_2005 <- .find_latest_grid_output("dis", "1691_1860", "1976_2005", all_months)
  gw_1976_2005 <- .find_latest_grid_output("rootmoist", "1691_1860", "1976_2005", all_months)
  case_basins <- read_sf("data/ardiv/selected_basins.gpkg")
  .run_plot02(bw_1976_2005, gw_1976_2005, "1976_2005", case_basins)

  bw_1931_1960 <- .find_latest_grid_output("dis", "1691_1860", "1931_1960", all_months)
  gw_1931_1960 <- .find_latest_grid_output("rootmoist", "1691_1860", "1931_1960", all_months)
  .run_plot02(bw_1931_1960, gw_1931_1960, "1931_1960")

}

