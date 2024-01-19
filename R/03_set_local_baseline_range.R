# ------------------------------------------------------------------------------
# 
# 03_set_local_baseline_range
# 
# Set the local baseline range (Fig. 1a) using baseline period quantiles given
# in config parameters local_bound_low and local_bound_high.
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

.get_quantile <- function(x, quant, cptmin, cptmax) {
  cpt <- x[3]
  cpt <- ifelse(cpt < cptmin | cpt > cptmax, cptmin, cpt)
  # x,y,cpt,quantile excluding cpt first data values
  return (c(x[1], x[2], quantile(x[(cpt+4):length(x)], quant)))
}

.determine_bounds <- function(file_in, cfg, baseline_ystart) {

  data <- readRDS(file_in)

  lblow <- data %>%
    select(x, y, cpt, starts_with(cfg$variable)) %>%
    as.matrix() %>%
    apply(MARGIN = 1, FUN = .get_quantile, quant = cfg$local_bound_low,
          cptmin = cfg$cptmin, cptmax = cfg$cptmax, simplify = FALSE) %>%
    bind_rows() %>%
    setNames(c("x", "y", "local_bound_low"))

  lbhigh <- data %>%
    select(x, y, cpt, starts_with(cfg$variable)) %>%
    as.matrix() %>%
    apply(MARGIN = 1, FUN = .get_quantile, quant = cfg$local_bound_high,
          cptmin = cfg$cptmin, cptmax = cfg$cptmax, simplify = FALSE) %>%
    bind_rows() %>%
    setNames(c("x", "y", "local_bound_high"))

  local_bounds <- data %>%
    select(-c(starts_with(cfg$variable), cpt, intercalibrated)) %>%
    left_join(lblow, by = c("x", "y")) %>%
    left_join(lbhigh, by = c("x", "y"))

  out <- file_in %>%
    str_replace_all(paste0("_monthly"), "_local_bounds")

  saveRDS(local_bounds, out)
  write(paste0("saved ", out), cfg$log_file, append = TRUE)

}

# public function
set_local_baseline_range <- function(cfg) {

  # setup
  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("creating dry and wet local bounds...", cfg$log_file, append = TRUE)
  vmon <- paste0(cfg$variable, "_monthly")
  ystart <- cfg$period_begin[which(cfg$period_label == cfg$baseline_period)]
  i_out <- cfg$impactmodel %>%
    str_replace_all("-", "_")
  f_out <- cfg$forcing %>%
    str_replace_all("-", "_")

  # check that enough directories exist
  local_bounds_fldr <- paste0("Data/", cfg$variable, "/", cfg$baseline_period,
                              "/", cfg$variable, "_local_bounds/")
  for (i in 1:length(i_out)) {
    for (j in 1:length(f_out)) {
      outdir <- paste0(local_bounds_fldr,
                       ifelse(i_out[i] == "", "", paste0(i_out[i], "/")),
                       f_out[j])
      if (!dir.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
    }
  }

  # query files to process
  files_in <- paste0("Data/", cfg$variable, "/",
                     cfg$baseline_period, "/", vmon) %>%
    list.files(full.names = TRUE, recursive = TRUE)
  i_out <- i_out %>% paste(collapse = "|")
  f_out <- f_out %>% paste(collapse = "|")
  files_in <- files_in[grepl(i_out, files_in) & grepl(f_out, files_in)]

  # compute and save local bounds
  for (i in 1:length(files_in)) {
    .determine_bounds(files_in[i], cfg, ystart)
  }
}
