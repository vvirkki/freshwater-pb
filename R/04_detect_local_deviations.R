# ------------------------------------------------------------------------------
# 
# 04_detect_local_deviations
# 
# Detect local deviations (Fig. 1b; also termed departures) in all grid cells,
# save deviations as differences to the local baseline range (dry deviations;
# values are < 0, wet deviations; values are > 0).
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

options(dplyr.summarise.inform = FALSE)

.check_departure <- function(obs, low, up) {

  d <- tibble(obs, low, up) %>%
    mutate(position = case_when(
      obs > up ~ obs - up,
      obs < low ~ obs - low
    ))
  return (d$position)

}

.count_dpts <- function(x, cptmin, cptmax, checklen, bound) {

  cpt <- x[3]
  cpt <- ifelse(cpt < cptmin | cpt > cptmax, cptmin, cpt)
  xvals <- x[(cpt+4):(cpt+4+checklen-1)]
  if (bound == "dry") {
    countdpts <- length(xvals[!is.na(xvals) & xvals < 0])
  } else if (bound == "wet") {
    countdpts <- length(xvals[!is.na(xvals) & xvals > 0])
  } else {
    stop("local bound not given correctly")
  }
  ret <- c(x[1], x[2], countdpts) %>%
    setNames(c("x", "y", paste0(bound, "Dpts")))
  return (ret)

}

.compare <- function(cfg, file_local_bounds, compare_to) {

  local_bounds <- readRDS(file_local_bounds)
  comparison_data <- file_local_bounds %>%
    str_replace_all("local_bounds", "monthly") %>%
    str_replace_all(cfg$baseline_period, compare_to) %>%
    readRDS() %>%
    as_tibble()

  data <- local_bounds %>%
    left_join(comparison_data, by = c("impactmodel", "forcing", "cell",
                                      "cellArea_km2", "x", "y"))

  departures <- data %>%
    mutate(across(starts_with(cfg$variable),
                  ~ .check_departure(.x, local_bound_low, local_bound_high))) %>%
    rename_with(~ paste0("dpt_", .x), starts_with(cfg$variable))

  missing_data_cells <- departures %>%
    filter(is.na(cpt)) %>%
    pull(cell)
  if (length(missing_data_cells) > 0) {
    write(paste0(length(missing_data_cells), " cells eliminated for no comparison data"),
          cfg$log_file, append = TRUE)
    departures <- departures %>%
      filter(!cell %in% missing_data_cells)
  }

  # save departures to tempdir if intercalibration is not complete yet
  if (compare_to == cfg$intcal_from & !cfg$is_corrected) {

    ic <- ifelse(cfg$intercalibrate, "corrected", "noncorrected")

    ref_dry_departures <- departures %>%
      select(x, y, cpt, starts_with("dpt_")) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .count_dpts, cptmin = cfg$cptmin, cptmax = cfg$cptmax,
            checklen = cfg$intcal_check_length, bound = "dry", simplify = FALSE) %>%
      bind_rows()

    ref_wet_departures <- departures %>%
      select(x, y, cpt, starts_with("dpt_")) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .count_dpts, cptmin = cfg$cptmin, cptmax = cfg$cptmax,
            checklen = cfg$intcal_check_length, bound = "wet", simplify = FALSE) %>%
      bind_rows()

    reference_departures <- ref_dry_departures %>%
      left_join(ref_wet_departures, by = c("x", "y")) %>%
      setNames(c("x", "y", paste0(ic, c("_refDryDpts", "_refWetDpts"))))

    out <- file_local_bounds %>%
      str_replace_all(cfg$baseline_period, "") %>%
      str_split("/") %>%
      unlist() %>%
      tail(1) %>%
      str_replace("local_bounds_", "ref_dpts") %>%
      paste0(cfg$tempdir, "/", ic, "_", .)

    saveRDS(reference_departures, out)
    write(paste0("saved ", out), cfg$log_file, append = TRUE)

  }

  return (departures)

}

# public function
detect_local_deviations <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("detecting local deviations...", cfg$log_file, append = TRUE)
  baseline_period <- cfg$baseline_period
  comparison_period <- cfg$period_label[which(cfg$period_label != baseline_period)]
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_")

  # check that enough directories exist
  for (p in 1:length(cfg$period_label)) {
    dpt_fldr <- paste0("Data/", cfg$variable, "/", cfg$period_label[p],
                       "/", cfg$variable, "_departures/")
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
  local_bounds <- paste0("Data/", cfg$variable, "/", baseline_period, "/",
                         cfg$variable, "_local_bounds") %>%
    list.files(recursive = TRUE, full.names = TRUE)
  i_out <- i_out %>% paste(collapse = "|")
  f_out <- f_out %>% paste(collapse = "|")
  local_bounds <- local_bounds[grepl(i_out, local_bounds) &
                               grepl(f_out, local_bounds)]

  if (!cfg$is_corrected & !cfg$intercalibrate) {
    write("saving deviations without intercalibration to tempdir...",
          cfg$log_file, append = TRUE)
  } else if (!cfg$is_corrected) {
    write("saving deviations after intercalibration to tempdir...",
          cfg$log_file, append = TRUE)
  }

  # compute and save departures
  for (i in 1:length(local_bounds)) {

    baseline_departures <- .compare(cfg, local_bounds[i], baseline_period)
    comparison_departures <- .compare(cfg, local_bounds[i], comparison_period)

    baseline_departures_out <- local_bounds[i] %>%
      str_replace_all("local_bounds", "departures")
    saveRDS(baseline_departures, baseline_departures_out)
    write(paste0("saved ", baseline_departures_out), cfg$log_file, append = TRUE)

    comparison_departures_out <- local_bounds[i] %>%
      str_replace_all("local_bounds", "departures") %>%
      str_replace_all(baseline_period, comparison_period)
    saveRDS(comparison_departures, comparison_departures_out)
    write(paste0("saved ", comparison_departures_out), cfg$log_file, append = TRUE)

  }
}
