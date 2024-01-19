# ------------------------------------------------------------------------------
# 
# output01_land_area_with_local_deviations
# 
# Determine pre-industrial variability (Fig. 1d) and create outputs (plot, csv)
# describing the percentage of land area with deviations, either globally or
# within areas prescribed in config.
#
# Figures using these outputs: 2, 4, 5
# Extended Data Figures using these outputs: 1, 2, 3, 4
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

.ts_boundVoc <- c("Summed dry and wet", "Dry", "Wet") %>%
  setNames(c("summed", "below", "above"))
.ts_varVoc <- c("streamflow", "root-zone soil moisture") %>%
  setNames(c("dis", "rootmoist"))

.read_visopts <- function(path) {

  # read config file; return a list of options to be used
  cfg_data <- readLines(path) %>%
    tibble() %>%
    rename(raw = 1) %>%
    separate(raw, into = c("variable", "value"), sep = ";")

  cfg <- cfg_data$value %>%
    as.list() %>%
    setNames(cfg_data$variable)

  for (i in 1:length(cfg)) {
    if (grepl(",", cfg[[i]])) {
      cfg[[i]] <- cfg[[i]] %>%
        str_split(",") %>%
        unlist() %>%
        as.vector()
    }
    suppressWarnings(
      cfg[[i]] <- ifelse(is.na(as.numeric(cfg[[i]])),
                         cfg[[i]], as.numeric(cfg[[i]]))
    )
  }

  return (cfg)

}

.compute_transgression_year <- function(boundary_values, area_shares,
                                        draw_class, npers) {

  for (i in 1:nrow(boundary_values)) {
    b <- boundary_values[i,]$boundary
    area_shares <- area_shares %>%
      mutate(!!b := boundary_values[i,]$value) %>%
      mutate(!!paste0("tg_", b) := landShareInClass > .[b])
  }

  cumulative_transgressions <- area_shares %>%
    mutate(across(starts_with("tg_"),
                  ~ zoo::rollsum(x = .x, k = npers, align = "right", fill = 0))) %>%
    arrange(year)

  persistent_transgression_years <- c() # year when transgression becomes persistent
  nm <- c()
  for (i in 1:nrow(boundary_values)) {
    b <- boundary_values[i,]$boundary
    tg_yr <- NA
    j <- 1
    while (j < nrow(cumulative_transgressions)) {
      yr <- cumulative_transgressions %>%
        slice(j:nrow(cumulative_transgressions)) %>%
        pull(!!paste0("tg_", b))
      if (all(yr == npers)) {
        tg_yr <- cumulative_transgressions %>%
          slice(j) %>%
          pull(year)
        persistent_transgression_years <- c(persistent_transgression_years, tg_yr)
        nm <- c(nm, paste0("tg_", b))
        break
      } else {
        j <- j + 1
      }
    }
    if (j == nrow(cumulative_transgressions)) {
      persistent_transgression_years <- c(persistent_transgression_years, tg_yr)
      nm <- c(nm, paste0("tg_", b))
    }
  }

  ret <- tibble(dptClass = draw_class,
                # first year of transgression that became persistent
                value = (persistent_transgression_years - npers + 1),
                boundary = str_replace(nm, "tg_", "persistent_transgression_year_")) %>%
    pivot_wider(id_cols = dptClass, names_from = boundary)

  return (ret)

}

.compose_output <- function(cfg, plt_data, draw_class, visopts) {

  # prepare data
  plt_data <- plt_data %>%
    filter(dptClass == draw_class)

  if (any(is.na(plt_data$landShareInClass))) {
    write("returning empty...", cfg$log_file, append = TRUE)
    return (list(ggplot(), tibble()) %>% setNames(c("plt", "currData")))
  }

  ystart <- max(min(cfg$period_begin) + cfg$cptmax, min(plt_data$year))
  yend <- max(plt_data$year)

  plt_data <- plt_data %>%
    filter(year >= ystart)
  period_duration <- plt_data %>%
    select(year, period) %>%
    distinct()

  # axis properties
  brky_max_fixed <- visopts[[paste0("ts_yrange_", draw_class)]]
  if (!as.logical(visopts$ts_fixed_yrange)) {
    brky_max <- max(brky_max_fixed, quantile(plt_data$landShareInClass, 0.95))
  } else {
    brky_max <- brky_max_fixed
  }
  breaks_y <- seq(0, brky_max, visopts$ts_yspacing) %>%
    round(2)

  breaks_x <- seq(ystart, yend, 1)
  breaks_x <- breaks_x[which(breaks_x %% visopts$ts_xspacing == 0)]
  breaks_x[1] <- ystart
  breaks_x[length(breaks_x)] <- yend

  # boundary lines and their labels
  boundaries <- plt_data %>%
    group_by(year) %>%
    summarise(landShareInClass = median(landShareInClass)) %>% # ensemble median
    ungroup() %>%
    left_join(period_duration, by = "year") %>%
    filter(period == cfg$baseline_period) %>%
    summarise(baseline = quantile(landShareInClass, cfg$boundary_baseline),
              upper_end = quantile(landShareInClass, cfg$boundary_upper_end)) %>%
    mutate(dptClass = draw_class)

  boundary_values <- data.frame(c("baseline", "upper_end"),
                                  c(cfg$boundary_baseline,
                                    cfg$boundary_upper_end)) %>%
    setNames(c("boundary", "position"))

  plt_data_boundaries <- boundaries %>%
    pivot_longer(-dptClass, names_to = "boundary") %>%
    mutate(xbegin = ystart,
           xend = yend)

  plt_data_boundary_labs <- plt_data_boundaries %>%
    left_join(boundary_values, by = "boundary") %>%
    mutate(label = paste0(sprintf("%.1f", round(value * 100, 1)), "%"),
           label = paste0(label, " (", position, ")"),
           xend = yend) %>%
    select(xend, value, label) %>%
    setNames(c("x", "y", "label"))

  # backgrounds
  bg_below_boundaries <- data.frame(space = "aaa_below_boundaries") %>%
    mutate(ymin = 0,
           ymax = min(plt_data_boundaries$value))

  bg_between_boundaries <- data.frame(space = "aab_between_boundaries") %>%
    mutate(ymin = min(plt_data_boundaries$value),
           ymax = max(plt_data_boundaries$value))

  bg_above_boundaries <- data.frame(space = "aac_above_boundaries") %>%
    mutate(ymin = max(plt_data_boundaries$value),
           ymax = Inf)

  plt_data_backgrounds <- bg_below_boundaries %>%
    bind_rows(bg_between_boundaries) %>%
    bind_rows(bg_above_boundaries) %>%
    mutate(xmin = ystart,
           xmax = yend)

  # interquartile range shadings
  plt_data_ribbons <- plt_data %>%
    group_by(year) %>%
    summarise(minrib = quantile(landShareInClass, 0.25),
              maxrib = quantile(landShareInClass, 0.75)) %>%
    ungroup() %>%
    left_join(period_duration, by = "year") %>%
    mutate(ribbonid = period)

  # lines and their labels
  plt_data_lines <- plt_data %>%
    group_by(year) %>%
    summarise(landShareInClass = median(landShareInClass)) %>% # ensemble median
    ungroup() %>%
    left_join(period_duration, by = "year") %>%
    mutate(lineid = paste0("median_", period))

  begin_rollmean <- max(cfg$period_begin) + cfg$cptmax - visopts$ts_moving_window_yrs + 1

  plt_data_lines_rollmean <- plt_data_lines %>%
    arrange(year) %>%
    filter(!(year >= max(cfg$period_begin) & year < begin_rollmean)) %>%
    group_by(lineid) %>%
    mutate(landShareInClass = c(rep(NA, visopts$ts_moving_window_yrs - 1),
                                zoo::rollmean(landShareInClass, visopts$ts_moving_window_yrs,
                                              align = "right"))) %>%
    filter(!is.na(landShareInClass)) %>%
    mutate(lineid = paste0("rollmean_", lineid)) %>%
    ungroup()

  plt_data_line_labs <- plt_data_lines_rollmean %>%
    filter(year == max(.$year)) %>%
    mutate(status = sprintf("%.1f", round(landShareInClass * 100, 1)),
           label = paste0(status, "% (", visopts$ts_moving_window_yrs, "-yr mean)")) %>%
    select(year, landShareInClass, label) %>%
    setNames(c("x", "y", "label"))

  dupl_row <- plt_data_lines %>%
    filter(year == max(cfg$period_begin) + cfg$cptmax) # to continue from dashed to solid line

  plt_data_lines <- plt_data_lines %>%
    mutate(dash = year >= max(cfg$period_begin) & year <= max(cfg$period_begin) + cfg$cptmax,
           lineid = ifelse(dash, paste0(lineid, "_dash"), lineid)) %>%
    select(-dash) %>%
    bind_rows(plt_data_lines_rollmean) %>%
    bind_rows(dupl_row)

  # annotations for persistent transgression times
  upper_end_boundary <- plt_data_boundaries %>%
    filter(boundary == "upper_end") %>%
    select(value, boundary)

  area_shares_rollmean <- plt_data_lines_rollmean %>%
    filter(period == cfg$period_label[cfg$period_label != cfg$baseline_period]) %>%
    select(year, landShareInClass)

  transgression_year <- .compute_transgression_year(upper_end_boundary,
                                                    area_shares_rollmean,
                                                    draw_class,
                                                    visopts$ts_persistent_transgression_yrs)

  plt_data_transgression_years <- transgression_year %>%
    pivot_longer(-dptClass, values_to = "year") %>%
    filter(!is.na(year)) %>%
    left_join(plt_data_lines_rollmean, by = "year") %>%
    select(year, landShareInClass) %>%
    rename(x = year, yend = landShareInClass) %>%
    mutate(xend = x, y = 0)

  # other plot elements
  plt_data_labs <- plt_data_boundary_labs %>%
    bind_rows(plt_data_line_labs)

  period_switch <- plt_data %>%
    group_by(period) %>%
    summarise(ymax = max(year)) %>%
    pull(ymax) %>%
    min()

  nensmem <- plt_data %>%
    group_by(impactmodel, forcing) %>%
    n_groups()

  pal_color <- c("black", "#545454", "black", "red", "red")
  pal_linetype <- c("solid", "dashed", "solid", "solid", "solid")
  plt_title <- paste0(.ts_boundVoc[[draw_class]], " ",
                      .ts_varVoc[[cfg$variable]],
                      " deviations (annual mean; ensemble n=",
                      nensmem, ")")

  # plot
  if (as.logical(visopts$ts_draw_backgrounds)) {
    plt <- ggplot() +
      geom_rect(data = plt_data_backgrounds,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = space),
                alpha = 0.3)
    pal_fill <- c("#3a4f8e", "#F4B326", "#cf5426", "grey50", "grey50")
  } else {
    plt <- ggplot()
    pal_fill <- c("grey50", "grey50")
  }

  plt <- plt +
    geom_ribbon(data = plt_data_ribbons,
                aes(x = year, ymin = minrib, ymax = maxrib, fill = ribbonid),
                alpha = 0.5) +
    geom_segment(data = plt_data_boundaries,
                 aes(x = xbegin, xend = xend, y = value, yend = value),
                 color = "grey25",
                 linetype = "dashed",
                 linewidth = (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_line(data = plt_data_lines,
              aes(x = year, y = landShareInClass, color = lineid, linetype = lineid),
              linewidth = 1.25 * (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_text(data = plt_data_transgression_years,
              aes(x = x, y = brky_max * 0.03, label = x),
              size = (visopts$ts_label_size / ggplot2::.pt),
              hjust = 0,
              nudge_x = 2) +
    geom_text(data = plt_data_labs,
              aes(x = x, y = y, label = label),
              size = (visopts$ts_label_size / ggplot2::.pt),
              hjust = 0) +
    geom_vline(xintercept = period_switch,
               color = "grey25",
               linewidth = (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_segment(aes(x = ystart, xend = yend, y = 0, yend = 0),
                 color = "grey25",
                 linewidth = (visopts$ts_guideline_width / ggplot2::.pt)) +
    geom_segment(data = plt_data_transgression_years,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "grey25",
                 linewidth = (visopts$ts_guideline_width / ggplot2::.pt) / 2) +
    scale_color_manual(values = pal_color) +
    scale_fill_manual(values = pal_fill) +
    scale_linetype_manual(values = pal_linetype) +
    scale_y_continuous(name = "Percentage of ice-free land area",
                       breaks = breaks_y,
                       labels = abs(breaks_y * 100),
                       expand = expansion()) +
    scale_x_continuous(breaks = breaks_x,
                       limits = c(ystart, yend),
                       expand = expansion(mult = c(0, 0.15))) +
    coord_cartesian(ylim = c(min(breaks_y), brky_max)) +
    ggtitle(plt_title) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = visopts$ts_title_size),
          axis.text = element_text(size = visopts$ts_axistext_size),
          axis.title = element_text(size = visopts$ts_annotation_size))

  # collect data on boundaries and current status for exporting
  current_status <- plt_data_lines_rollmean %>%
    filter(year == yend) %>%
    select(landShareInClass) %>%
    mutate(boundary = "current") %>%
    rename(value = landShareInClass)

  export_boundaries_current_status <- plt_data_boundaries %>%
    filter(boundary %in% c("baseline", "upper_end")) %>%
    select(value, boundary) %>%
    bind_rows(current_status) %>%
    mutate(dptClass = draw_class) %>%
    pivot_wider(id_cols = dptClass, names_from = boundary) %>%
    mutate(relStatus_baseline = current / baseline,
           relStatus_upper_end = current / upper_end,
           area = unique(plt_data$area)) %>%
    left_join(transgression_year, by = "dptClass") %>%
    rename(class = dptClass) %>%
    mutate(class = case_when(
             draw_class == "summed" ~ "summed_deviations",
             draw_class == "above" ~ "wet_deviations",
             draw_class == "below" ~ "dry_deviations")) %>%
    select(area, everything())

  export_ensemble_median_iqr <- plt_data_lines %>%
    filter(!grepl("rollmean", lineid)) %>%
    filter(!(lineid == "median_historical_dash" & year == 1891)) %>%
    select(year, landShareInClass) %>%
    mutate(class = case_when(
             draw_class == "summed" ~ "summed_deviations",
             draw_class == "above" ~ "wet_deviations",
             draw_class == "below" ~ "dry_deviations"
           ),
           area = unique(plt_data$area)) %>%
    left_join(plt_data_ribbons %>% select(year, minrib, maxrib), by = "year") %>%
    rename(ensemble_median = landShareInClass,
           IQR_min = minrib,
           IQR_max = maxrib) %>%
    select(area, year, class, everything()) %>%
    arrange(year)

  return (list(plt,
               export_boundaries_current_status,
               export_ensemble_median_iqr) %>%
            setNames(c("plt",
                       "export_boundaries_current_status",
                       "export_ensemble_median_iqr")))

}

# public function
output_land_area_with_local_deviations <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("outputting land area with local deviations...",
        cfg$log_file, append = TRUE)

  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")
  ardiv_lbl <- ifelse(cfg$areal_division == "", "global", cfg$areal_division)

  visopts_cfg <- paste0("configs/", cfg$visopts_cfg)
  if (!file.exists(visopts_cfg)) {
    stop("visopts not defined for this config")
  } else {
    visopts <- .read_visopts(visopts_cfg)
  }

  if (!cfg$areal_division == "") {
    unique_areas <- paste0("data/ardiv/", cfg$areal_division, ".tif") %>%
      rast() %>%
      values() %>%
      unique() %>%
      as.numeric()
    unique_areas <- unique_areas[!is.na(unique_areas)]
  } else {
    unique_areas <- c(NA)
  }

  # create directories
  tstamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_fldr <- paste0("output/", cfg$variable, "_", ardiv_lbl,
                     "_land_area_with_local_deviations_",
                     paste(cfg$impactmodel, collapse = "_"), "_", tstamp)
  if (!dir.exists(output_fldr)) {
    dir.create(output_fldr, recursive = TRUE)
  }
  if (any(!is.na(unique_areas))) {
    for (i in 1:length(unique_areas)) {
      figs_area_out <- paste0(output_fldr, "/id", unique_areas[i])
      if (!dir.exists(figs_area_out)) {
        dir.create(figs_area_out, recursive = TRUE)
      }
    }
  }

  # query files to process
  aggr_files <- list.files(paste0("Data/", cfg$variable),
                           recursive = TRUE, full.names = TRUE)
  aggr_files <- aggr_files[grepl(paste0(ardiv_lbl, "_departure_aggregates"), aggr_files) &
                           grepl(f_out, aggr_files) &
                           grepl(i_out, aggr_files)]
  aggr_data <- aggr_files %>%
    lapply(readRDS) %>%
    bind_rows()

  # unique ensemble members in all data
  nensmem <- aggr_data %>%
    group_by(impactmodel, forcing) %>%
    n_groups()

  # create output for each area
  for (i in 1:length(unique_areas)) {

    curr_status_data <- tibble()
    if (!is.na(unique_areas[i])) {
      write(paste0("processing area id ", unique_areas[i], "..."),
            cfg$log_file, append = TRUE)
      area_id <- paste0("id", unique_areas[i])
      filter_area <- paste0(cfg$areal_division, "_", area_id)
    } else {
      area_id <- "global"
      filter_area <- "global"
      write("processing global...", cfg$log_file, append = TRUE)
    }

    area_data <- aggr_data %>%
      filter(area == filter_area)
    area_fldr <- ifelse(!is.na(unique_areas[i]), paste0("id", unique_areas[i], "/"), "")

    # check if an area is covered by all or not all ensemble members
    ngroups <- area_data %>%
      group_by(impactmodel, forcing) %>%
      n_groups()

    if (ngroups == 0) {
      write("some small area, all ensemble members missing...",
            cfg$log_file, append = TRUE)
      next
    } else if (ngroups < nensmem) {
      write(paste0(ngroups, " ensemble members covering area..."),
            cfg$log_file, append = TRUE)
    }

    # take annual means of months
    aggregates_year <- area_data %>%
      filter(!is.na(dptClass)) %>%
      group_by(variable, impactmodel, forcing, area, period, year, dptClass) %>%
      summarise(landShareInClass = mean(landShareInClass)) %>%
      ungroup()

    aggregates_year_summed <- aggregates_year %>%
      filter(!is.na(dptClass)) %>%
      group_by(variable, impactmodel, forcing, area, period, year) %>%
      summarise(landShareInClass = sum(landShareInClass)) %>%
      ungroup() %>%
      mutate(dptClass = "summed")

    output_dry_dpts <- .compose_output(cfg, aggregates_year, "below", visopts)
    output_wet_dpts <- .compose_output(cfg, aggregates_year, "above", visopts)
    output_summed_dpts <- .compose_output(cfg, aggregates_year_summed, "summed", visopts)

    arr_annual <- gridExtra::arrangeGrob(grobs = list(output_summed_dpts$plt,
                                                      output_wet_dpts$plt,
                                                      output_dry_dpts$plt))

    plt_year_out <- paste0(output_fldr, "/", area_fldr, cfg$variable, "_",
                           ifelse(!cfg$areal_division == "", paste0(cfg$areal_division, "_"), ""),
                                  area_id,
                           "_land_area_with_local_deviations_annual_mean.pdf")
    ggsave(plt_year_out, arr_annual, width = visopts$ts_out_width,
           height = 3 * visopts$ts_out_height, units = visopts$ts_out_units)
    write(paste0("saved ", plt_year_out, "..."), cfg$log_file, append = TRUE)

    export_boundaries <- output_summed_dpts$export_boundaries_current_status %>%
      bind_rows(output_wet_dpts$export_boundaries_current_status) %>%
      bind_rows(output_dry_dpts$export_boundaries_current_status)

    export_lines <- output_summed_dpts$export_ensemble_median_iqr %>%
      bind_rows(output_wet_dpts$export_ensemble_median_iqr) %>%
      bind_rows(output_dry_dpts$export_ensemble_median_iqr)

    export_boundaries_out <- paste0(output_fldr, "/", area_fldr, cfg$variable, "_",
                                    ifelse(!cfg$areal_division == "", paste0(cfg$areal_division, "_"), ""),
                                    area_id,
                                    "_boundaries_current_status_persistent_transgression_year.csv")
    write.csv(export_boundaries, export_boundaries_out, row.names = FALSE)
    write(paste0("saved ", export_boundaries_out, "..."), cfg$log_file, append = TRUE)

    export_lines_out <- paste0(output_fldr, "/", area_fldr, cfg$variable, "_",
                               ifelse(!cfg$areal_division == "", paste0(cfg$areal_division, "_"), ""),
                               area_id,
                               "_land_area_with_local_deviations_annual_mean_ensemble_median_IQR.csv")
    write.csv(export_lines, export_lines_out, row.names = FALSE)
    write(paste0("saved ", export_lines_out, "..."), cfg$log_file, append = TRUE)

  }
}
