# ------------------------------------------------------------------------------
# 
# main
# 
# Main analysis script that calls for required libraries, initialises analysis
# and output functions, prepares and calls configs for parallel data
# pre-processing and creates outputs used in the manuscript.
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

# Load libraries & scripts -----------------------------------------------------

library(qmap)
library(terra)
library(sf)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)
library(tmap)
library(viridisLite)
library(changepoint)
library(zoo)
library(doParallel)

# data preparation scripts
source("R/00_configure_run.R")
source("R/01_prepare_monthly_data.R")
source("R/02_intercalibrate_periods.R")
source("R/03_set_local_baseline_range.R")
source("R/04_detect_local_deviations.R")
source("R/05_build_ice_areas.R")
source("R/06_get_land_area_with_local_deviations.R")
source("R/07_get_local_deviation_frequency.R")

# data export & plotting scripts
source("R/output01_land_area_with_local_deviations.R")
source("R/output02_local_deviation_frequency_change.R")
source("R/plot01_regional_persistent_transgressions.R")
source("R/plot02_combined_local_deviation_frequency_increases.R")

# save sessionInfo()
sink(paste0("sessionInfo ", Sys.time(), ".txt"))
sessionInfo()
sf::sf_extSoftVersion()
sink()

# Prepare configs for parallel running -----------------------------------------

# template config
cfg <- configure_run("configs/template main isimip2")

# variables
variables <- c("rootmoist", "dis")

# ISIMIP2 GCMs
gcms <- c("gfdl-esm2m", "hadgem2-es", "ipsl-cm5a-lr", "miroc5")

parallel_configs <- list()

for (v in 1:length(variables)) {

  if (variables[v] == "dis") {
    ghms <- c("h08", "lpjml", "pcr-globwb", "watergap2", "mpi-hm", "matsiro")
  } else if (variables[v] == "rootmoist") {
    ghms <- c("clm50", "lpjml", "mpi-hm", "pcr-globwb")
  }

  for (i in 1:length(ghms)) {
    for (j in 1:length(gcms)) {

      ensemble_member <- paste0(variables[v], " ", ghms[i], " ", gcms[j])
      new_config <- cfg
      new_config$variable <- variables[v]
      new_config$impactmodel <- ghms[i]
      new_config$forcing <- gcms[j]

      if ((ghms[i] == "mpi-hm") & gcms[j] == "hadgem2-es") {
        # missing ensemble member
        next
      } else {
        # open log file
        if (!dir.exists("logs")) { dir.create("logs") }
        new_config$log_file <- paste0("logs/log ", Sys.time(),
                                      " ", ensemble_member, ".txt")
        write(format(Sys.time(), "%a %b %d %X %Y"), new_config$log_file, append = TRUE)
        write("opening log...", new_config$log_file, append = TRUE)

        parallel_configs[[length(parallel_configs) + 1]] <- new_config
        message(paste0("prepared config: ", ensemble_member))

      }
    }
  }
}

# avoid messing up with globalenv in parallel processing
unlink(cfg$log_file)
remove(cfg, ensemble_member, gcms, ghms, i, j, new_config)

# Prepare main analysis data with prescribed configs ---------------------------

registerDoParallel(cores = detectCores() / 2)
foreach(i = 1:length(parallel_configs),
        .errorhandling = "pass",
        .verbose = TRUE) %dopar% {

          # prepare data
          prepare_monthly_data(parallel_configs[[i]])

          # run until detecting local deviations without intercalibration
          parallel_configs[[i]]$intercalibrate <- FALSE
          intercalibrate_periods(parallel_configs[[i]])
          set_local_baseline_range(parallel_configs[[i]])
          detect_local_deviations(parallel_configs[[i]])

          # do intercalibration for all cells and run until local deviations
          parallel_configs[[i]]$intercalibrate <- TRUE
          intercalibrate_periods(parallel_configs[[i]])
          set_local_baseline_range(parallel_configs[[i]])
          detect_local_deviations(parallel_configs[[i]])

          # check if intercalibration has improved the fit between periods
          parallel_configs[[i]]$is_corrected <- TRUE
          intercalibrate_periods(parallel_configs[[i]])

          # create the final local deviations using the best possible intercalibration
          set_local_baseline_range(parallel_configs[[i]])
          detect_local_deviations(parallel_configs[[i]])

          # build ice areas (cells that are permanent ice)
          build_ice_areas(parallel_configs[[i]])

          # global percentages of land areas with local deviations
          get_land_area_with_local_deviations(parallel_configs[[i]])

          # local deviation frequencies for 1691-1860 and
          # 1976-2005 (default config)
          get_local_deviation_frequency(parallel_configs[[i]])

          # 1931-1960
          parallel_configs[[i]]$localfreq_period_begin <- c(1691, 1931)
          parallel_configs[[i]]$localfreq_period_end <- c(1860, 1960)
          get_local_deviation_frequency(parallel_configs[[i]])

          # hybas2 regional percentages of land areas with local deviations
          parallel_configs[[i]]$areal_division <- "hybas2"
          get_land_area_with_local_deviations(parallel_configs[[i]])

          # hybas3 regional percentages of land areas with local deviations
          parallel_configs[[i]]$areal_division <- "hybas3"
          get_land_area_with_local_deviations(parallel_configs[[i]])

        }
stopImplicitCluster()

# Create streamflow outputs ----------------------------------------------------

cfg <- configure_run("configs/dis main isimip2")

# Global land area with local deviations (Fig. 2, Extended Data Fig. 1)
output_land_area_with_local_deviations(cfg)

# GHM medians of global land area with local deviations (Extended Data Fig. 3)
ghms <- c("h08", "lpjml", "pcr-globwb", "watergap2", "mpi-hm", "matsiro")
for (i in 1:length(ghms)) {
  cfg$impactmodel <- ghms[i]
  output_land_area_with_local_deviations(cfg)
}
cfg$impactmodel <- ghms # set all GHMs back to config

# Regional land area with local deviations (Fig. 4)
cfg$areal_division <- "hybas2"
output_land_area_with_local_deviations(cfg)
cfg$areal_division <- "hybas3"
output_land_area_with_local_deviations(cfg)

# Local deviation frequency (Fig. 3, 5, Extended Data Fig. 5)
output_local_deviation_frequency_change(cfg)

# Seasonal local deviation frequency (Extended Data Fig. 6)
cfg$localfreq_season <- c("December", "January", "February")
output_local_deviation_frequency_change(cfg)
cfg$localfreq_season <- c("June", "July", "August")
output_local_deviation_frequency_change(cfg)

# Local deviation frequency in 1931-1960 (Extended Data Fig. 7)
cfg$localfreq_season <- ""
cfg$localfreq_period_begin <- c(1691, 1931)
cfg$localfreq_period_end <- c(1860, 1960)
output_local_deviation_frequency_change(cfg)

# Create soil moisture outputs -------------------------------------------------

cfg <- configure_run("configs/rootmoist main isimip2")

# Global land area with local deviations (Fig. 2, Extended Data Fig. 1)
output_land_area_with_local_deviations(cfg)

# GHM medians of global land area with local deviations (Extended Data Fig. 4)
ghms <- c("clm50", "lpjml", "mpi-hm", "pcr-globwb")
for (i in 1:length(ghms)) {
  cfg$impactmodel <- ghms[i]
  output_land_area_with_local_deviations(cfg)
}
cfg$impactmodel <- ghms # set all GHMs back to config

# Regional land area with local deviations (Fig. 4)
cfg$areal_division <- "hybas2"
output_land_area_with_local_deviations(cfg)
cfg$areal_division <- "hybas3"
output_land_area_with_local_deviations(cfg)

# Local deviation frequency (Fig. 3, 5, Extended Data Fig. 5)
output_local_deviation_frequency_change(cfg)

# Seasonal local deviation frequency (Extended Data Fig. 6)
cfg$localfreq_season <- c("December", "January", "February")
output_local_deviation_frequency_change(cfg)
cfg$localfreq_season <- c("June", "July", "August")
output_local_deviation_frequency_change(cfg)

# Local deviation frequency in 1931-1960 (Extended Data Fig. 7)
cfg$localfreq_season <- ""
cfg$localfreq_period_begin <- c(1691, 1931)
cfg$localfreq_period_end <- c(1860, 1960)
output_local_deviation_frequency_change(cfg)

# Create combined streamflow and soil moisture outputs -------------------------

# put Fig. 4 together; not run with config, above data preparation required
plot_regional_persistent_transgressions()

# put Fig. 5 together; not run with config, above data preparation required
plot_combined_local_deviation_frequency_increases()

# Detected changepoints and cells with accepted intercalibration ---------------

.check_cpt <- function(x) {
  d <- readRDS(x)
  return (nrow(d %>% filter(cpt > 10 & cpt <= 30)) / nrow(d))
}

.check_intercalibration <- function(x) {
  d <- readRDS(x)
  return (nrow(d %>% filter(intercalibrated)) / nrow(d))
}

fldrs <- list("data/dis/preind/dis_monthly",
              "data/dis/historical/dis_monthly",
              "data/rootmoist/preind/rootmoist_monthly",
              "data/rootmoist/historical/rootmoist_monthly") %>%
  setNames(c("dis_preind", "dis_hist", "rootmoist_preind", "rootmoist_hist"))

# iterate through monthly data files by all ensemble members
cpt_detections <- fldrs %>%
  lapply(list.files, recursive = TRUE, full.names = TRUE) %>%
  lapply(lapply, .check_cpt) %>%
  lapply(unlist)

cpt_detections %>% lapply(min)
cpt_detections %>% lapply(max)

accepted_intercalibrations <- fldrs %>%
  lapply(list.files, recursive = TRUE, full.names = TRUE) %>%
  lapply(lapply, .check_intercalibration) %>%
  lapply(unlist)

accepted_intercalibrations %>% lapply(min)
accepted_intercalibrations %>% lapply(max)
