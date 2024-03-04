# ------------------------------------------------------------------------------
# 
# 02_intercalibrate_periods
# 
# Perform iterative intercalibration between baseline and comparison periods
# described in config; example shown in Extended Data Fig. 9.
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

.extrapolate_linear <- function(x, cptmin, cptmax) {

  # setup
  cpt <- x[3]
  cpt <- ifelse(cpt < cptmin | cpt > cptmax, cptmin, cpt)
  fit_data <- x[(4+cpt):length(x)]
  lf <- length(fit_data)

  # model + normal noise
  lmodel <- lm(y ~ x, data.frame(1:lf, fit_data) %>% setNames(c("x", "y")))
  normal_noise <- rnorm(100, mean = 0, sd = sd(fit_data))

  # extrapolation
  extrapolate_data <- data.frame((lf+1):(lf+100)) %>% setNames("x")
  extrapolation <- predict(lmodel, extrapolate_data) + normal_noise

  return (c(x[1], x[2], extrapolation))

}

.correct_qmap <- function(x, obslen, modlen, corrlen, trlen, nq = 4) {

  trbegin <- x[3]+1
  excl_corr <- x[4]

  obs_data_full <- x[5:(5+obslen-1)]
  mod_data_full <- x[(5+obslen):(5+obslen+modlen-1)]
  tocorr_data_excluded <- x[(5+obslen+modlen):(5+excl_corr+obslen+modlen-1)]
  tocorr_data <- x[(5+excl_corr+obslen+modlen):(5+obslen+modlen+corrlen-1)]
  nm_tocorr <- c(names(tocorr_data_excluded), names(tocorr_data))

  obs_data <- obs_data_full[trbegin:(trbegin+trlen-1)]
  mod_data <- mod_data_full[trbegin:(trbegin+trlen-1)]

  ret <- tryCatch({
    # attempt fitting and applying qmap
    qmap_fit <- fitQmapRQUANT(obs = obs_data,
                              mod = mod_data,
                              wet.day = FALSE,
                              qstep = 1 / nq)
    qmap_apply <- doQmapRQUANT(tocorr_data, qmap_fit)
    ret <- c(FALSE, tocorr_data_excluded, qmap_apply) %>%
      setNames(c("warn", nm_tocorr))
  }, warning = function(cond) {
    if (any(obs_data > 0) & any(mod_data > 0)) {
      if (var(obs_data) > 0 & var(mod_data) > 0) {
        # qmap outputs a warning with very small values --> mark as failing
        ret <- c(TRUE, rep(-1, corrlen)) %>%
          setNames(c("warn", nm_tocorr))
      } else {
        # all values are non-zero but constant with no variance --> qmap fails
        ret <- c(FALSE, rep(-1, corrlen)) %>%
          setNames(c("warn", nm_tocorr))
      }
    } else {
      # all values are zero --> qmap fails
      ret <- c(FALSE, rep(-1, corrlen)) %>%
        setNames(c("warn", nm_tocorr))
    }
  })
  return(c(x[1], x[2], ret))

}

.do_intercalibration <- function(file_in, cfg) {

  # intercalibration begins
  data_to_correct <- readRDS(file_in) %>%
    as_tibble() %>%
    rename(to_cpt = cpt)

  # no intercalibration done --> perform correction for all data
  if (!cfg$is_corrected) {

    data_correction_from <- file_in %>%
      str_replace_all(cfg$intcal_to, cfg$intcal_from) %>%
      readRDS() %>%
      as_tibble() %>%
      rename(from_cpt = cpt)

    # linear extrapolation
    extrapolation_100yrs <- data_to_correct %>%
      select(x, y, to_cpt, starts_with(cfg$variable)) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .extrapolate_linear,
            cptmin = cfg$cptmin, cptmax = cfg$cptmax, simplify = FALSE)%>%
      bind_rows() %>%
      setNames(c("x", "y", paste0("expol", 1:100)))

    data_to_correct_extrapolated <- data_to_correct %>%
      left_join(extrapolation_100yrs, by = c("x" = "x", "y" = "y"))

    # quantile mapping correction
    qmap_data_to_correct <- data_to_correct %>%
      mutate(to_fci = ifelse(to_cpt < cfg$cptmin | to_cpt > cfg$cptmax,
                             cfg$cptmin, to_cpt)) %>%
      select(x, y, to_fci, starts_with(cfg$variable))
    ndata_to_correct <- ncol(qmap_data_to_correct) - 3

    qmap_fit_data_mod <- data_to_correct_extrapolated %>%
      select(x, y, starts_with("expol"))
    ndata_mod <- ncol(qmap_fit_data_mod) - 2

    qmap_fit_data_obs <- data_correction_from %>%
      mutate(from_fci = ifelse(from_cpt < cfg$cptmin | from_cpt > cfg$cptmax,
                               cfg$cptmin, from_cpt)) %>%
      select(x, y, from_fci, starts_with(cfg$variable))
    ndata_obs <- ncol(qmap_fit_data_obs) - 3

    qmap_correction_data <- qmap_fit_data_obs %>%
      left_join(qmap_fit_data_mod, by = c("x", "y")) %>%
      left_join(qmap_data_to_correct, by = c("x", "y")) %>%
      select(x, y, from_fci, to_fci, everything()) %>%
      as.matrix() %>%
      apply(MARGIN = 1, FUN = .correct_qmap,
            obslen = ndata_obs, modlen = ndata_mod, corrlen = ndata_to_correct,
            trlen = cfg$intcal_period_length, nq = cfg$qmap_nquantiles,
            simplify = FALSE) %>%
      bind_rows()

    qmap_corrected_data <- data_to_correct %>%
      select(-starts_with(cfg$variable)) %>%
      left_join(qmap_correction_data, by = c("x", "y"))

    qmap_failing <- qmap_corrected_data %>%
      filter(if_any(starts_with(cfg$variable), ~ .x == -1))

    if (nrow(qmap_failing) > 0) {
      # log the amount of cells in which qmap fails and warns
      write(paste0("qmap failed in ", nrow(qmap_failing), " instances, ",
                   qmap_failing %>% filter(warn == 1) %>% nrow(),
                   " of which are due to warning"),
            cfg$log_file, append = TRUE)
      qmap_failing_cells <- qmap_failing$cell
    } else {
      qmap_failing_cells <- c()
    }

    uncorrected_data <- data_to_correct %>%
      filter(cell %in% qmap_failing_cells)

    corrected_data <- qmap_corrected_data %>%
      select(-warn) %>%
      filter(!cell %in% qmap_failing_cells) %>%
      mutate(intercalibrated = TRUE)

    data_out <- corrected_data %>%
      bind_rows(uncorrected_data) %>%
      arrange(cell) %>%
      rename(cpt = to_cpt)

  } else { # all data has been intercalibrated --> check if departures decrease

    refs_from <- file_in %>%
      str_split("/") %>%
      unlist() %>%
      tail(1)

    ref_dpts_nocorr <- refs_from %>%
      str_replace(paste0(cfg$variable, "_monthly_", cfg$intcal_to),
                  paste0("noncorrected_", cfg$variable, "_ref_dpts")) %>%
      paste0(cfg$tempdir, "/", .) %>%
      readRDS()

    ref_dpts_corr <- refs_from %>%
      str_replace(paste0(cfg$variable, "_monthly_", cfg$intcal_to),
                  paste0("corrected_", cfg$variable, "_ref_dpts")) %>%
      paste0(cfg$tempdir, "/", .) %>%
      readRDS()

    check_ref_dpts <- ref_dpts_nocorr %>%
      left_join(ref_dpts_corr, by = c("x", "y")) %>%
      mutate(change_wetDpts = corrected_refWetDpts - noncorrected_refWetDpts,
             change_dryDpts = corrected_refDryDpts - noncorrected_refDryDpts,
             change_net = change_wetDpts + change_dryDpts)

    correction_result <- check_ref_dpts %>%
      mutate(correction_failing = ifelse(change_net >= 0, TRUE, FALSE)) %>%
      select(x, y, correction_failing)

    uncorrected_data <- refs_from %>%
      paste0(cfg$tempdir, "/", .) %>%
      readRDS() %>%
      left_join(correction_result, by = c("x", "y")) %>%
      filter(correction_failing) %>%
      select(-correction_failing)

    # do some sanity checks
    if (any(uncorrected_data$intercalibrated)) {
      stop("data intercalibrated (shouldn't be)")
    }

    corrected_data <- data_to_correct %>%
      left_join(correction_result, by = c("x", "y")) %>%
      filter(!correction_failing) %>%
      rename(cpt = to_cpt) %>%
      select(-correction_failing)

    if (any(!corrected_data$intercalibrated)) {
      stop("data not intercalibrated (should be)")
    }

    data_out <- corrected_data %>%
      bind_rows(uncorrected_data) %>%
      arrange(cell)

    if (nrow(uncorrected_data) + nrow(corrected_data) != nrow(data_to_correct)) {
      stop("mismatch in row numbers")
    }

    # write log outputs
    write(paste0("saving ", data_out %>% filter(!intercalibrated) %>% nrow(), " cells uncorrected..."),
          cfg$log_file, append = TRUE)
    write(paste0("saving ", data_out %>% filter(intercalibrated) %>% nrow(), " cells corrected..."),
          cfg$log_file, append = TRUE)

  }

  saveRDS(data_out, file_in) # replace old data
  write(paste0("saved ", file_in), cfg$log_file, append = TRUE)

}

# public function
intercalibrate_periods <- function(cfg) {

  # setup
  set.seed(42)

  vmon <- paste0(cfg$variable, "_monthly")
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_") %>%
    paste(collapse = "|")

  # query files to process
  files_in <- paste0("Data/", cfg$variable, "/", cfg$intcal_to, "/", vmon) %>%
    list.files(full.names = TRUE, recursive = TRUE)
  files_in <- files_in[grepl(i_out, files_in) & grepl(f_out, files_in)]

  # when not doing intercalibration, save uncorrected data to tempdir
  if (!as.logical(cfg$intercalibrate)) {
    write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
    write("no intercalibration requested, saving monthly data to tempdir...",
          cfg$log_file, append = TRUE)
    for (i in 1:length(files_in)) {
      tempfile_out <- files_in[i] %>%
        str_split("/") %>%
        unlist() %>%
        tail(1) %>%
        paste0(cfg$tempdir, "/", .)
      file.copy(files_in[i], tempfile_out)
      write(paste0("saved ", tempfile_out), cfg$log_file, append = TRUE)
    }
  } else {
    if (!cfg$is_corrected) {
      write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
      write("performing intercalibration for all data...",
            cfg$log_file, append = TRUE)
    } else {
      write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
      write("checking if intercalibration improves fit between periods...",
            cfg$log_file, append = TRUE)
    }
    # compute and save bias corrected data
    for (i in 1:length(files_in)) {
      .do_intercalibration(files_in[i], cfg)
    }
  }

}
