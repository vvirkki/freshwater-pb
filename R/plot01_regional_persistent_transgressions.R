# ------------------------------------------------------------------------------
# 
# plot01_regional_persistent_transgressions
# 
# Find years of persistent transgression of the upper end of pre-industrial
# variability in areas prescribed in config and plot vector maps of the
# transgression year.
#
# Figures using these outputs: 4
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

# the latest output folder with given variable and areal division
.find_latest_regional_output <- function(variable_ardiv) {

  opf <- list.files("output", full.names = TRUE)

  latest_fldr <- opf[grepl(variable_ardiv, opf)] %>%
    tibble(fldr = .) %>%
    mutate(ctime = file.info(fldr)$ctime) %>%
    arrange(desc(ctime)) %>%
    slice(1) %>%
    pull(fldr)

  return (latest_fldr)

}

.draw_map <- function(geoms, fldr, draw_variable, sum_only = TRUE) {

  files_in <- list.files(fldr, recursive = TRUE, full.names = TRUE)
  files_in <- files_in[grepl("boundaries_current_status_persistent_transgression_year", files_in)]

  data_in <- files_in %>%
    lapply(read_csv, show_col_types = FALSE) %>%
    bind_rows() %>%
    rename(transgression_year = !!draw_variable) %>%
    mutate(transgression_year = transgression_year + 9) # to comply with Fig. 4 definition

  plt_data <- data_in %>%
    left_join(geoms, by = c("area" = "id")) %>%
    st_as_sf() %>%
    tidyr::separate(area, into = c("hybas", "id")) %>%
    mutate(id = readr::parse_number(id),
           idfirst = substr(as.character(id), 1, 1),
           transgression_year = ifelse(as.numeric(idfirst) == 9, -1, transgression_year)) %>% # hybases in Greenland
    select(hybas, id, class, transgression_year)

  missing_geoms <- geoms %>%
    filter(!id %in% data_in$area) %>%
    slice(rep(1:n(), 3)) %>%
    arrange(id) %>%
    tidyr::separate(id, into = c("hybas", "id")) %>%
    mutate(class = rep(c("summed_deviations", "wet_deviations", "dry_deviations"),
                       nrow(.) / 3),
           transgression_year = NA,
           id = readr::parse_number(id)) %>%
    select(hybas, id, class, transgression_year)

  plt_data <- plt_data %>%
    bind_rows(missing_geoms)

  if (sum_only) {
    plt_data <- plt_data %>%
      filter(class == "summed_deviations")
  }

  brk <- c(-Inf, 0, 1940, 1950, 1960, 1970, 1980, 1990, Inf)
  pal <- c("#b3b3b3", viridisLite::plasma(7, begin = 0.2, end = 0.8, direction = -1))
  lbl <- paste0(brk[-length(brk)], "-", brk[-1])
  lbl[1] <- "ice"
  lbl[2] <- str_replace(lbl[2], "0-", "-")
  lbl[length(lbl)] <- str_replace(lbl[length(lbl)], "-Inf", "-")

  mm <- tm_shape(plt_data, projection = "ESRI:54030") +
    tm_polygons(col = "transgression_year",
                title = "persistent transgression year",
                breaks = brk,
                palette = pal,
                labels = lbl,
                midpoint = 0,
                textNA = "no transgression",
                colorNA = "#e5e5e5",
                lwd = 0.1,
                legend.is.portrait = FALSE) +
    tm_facets(by = "class") +
    tm_layout(frame = FALSE,
              legend.outside = TRUE,
              legend.outside.position = "bottom",
              legend.text.size = 0.2,
              legend.title.size = 0.2)

  return (mm)

}

.run_plot01 <- function() {

    if (!dir.exists("output/regional_persistent_transgression_plots")) {
      dir.create("output/regional_persistent_transgression_plots", recursive = TRUE)
    }

    # geometries
    hybas2 <- read_sf("data/ardiv/hybas2.gpkg") %>%
      mutate(id = paste0("hybas2_id", id))
    hybas3 <- read_sf("data/ardiv/hybas3.gpkg") %>%
      mutate(id = paste0("hybas3_id", id))

    # aggregate plots & saved transgression time data folders
    bw_hybas2 <- .find_latest_regional_output("dis_hybas2")
    bw_hybas3 <- .find_latest_regional_output("dis_hybas3")
    gw_hybas2 <- .find_latest_regional_output("rootmoist_hybas2")
    gw_hybas3 <- .find_latest_regional_output("rootmoist_hybas3")

    w <- 70
    h <- 40

    map_bw_hybas2 <- .draw_map(hybas2, bw_hybas2,  "persistent_transgression_year_upper_end")
    tmap_save(map_bw_hybas2, "output/regional_persistent_transgression_plots/dis_persistent_transgression_year_hybas2.pdf",
              width = w, height = h, units = "mm")

    map_bw_hybas3 <- .draw_map(hybas3, bw_hybas3, "persistent_transgression_year_upper_end")
    tmap_save(map_bw_hybas3, "output/regional_persistent_transgression_plots/dis_persistent_transgression_year_hybas3.pdf",
              width = w, height = h, units = "mm")

    map_gw_hybas2 <- .draw_map(hybas2, gw_hybas2, "persistent_transgression_year_upper_end")
    tmap_save(map_gw_hybas2, "output/regional_persistent_transgression_plots/rootmoist_persistent_transgression_year_hybas2.pdf",
              width = w, height = h, units = "mm")

    map_gw_hybas3 <- .draw_map(hybas3, gw_hybas3, "persistent_transgression_year_upper_end")
    tmap_save(map_gw_hybas3, "output/regional_persistent_transgression_plots/rootmoist_persistent_transgression_year_hybas3.pdf",
              width = w, height = h, units = "mm")

}

# public function
plot_regional_persistent_transgressions <- function() {
  .run_plot01()
}
