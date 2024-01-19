# ------------------------------------------------------------------------------
# 
# sensitivity_analysis
# 
# Script in the style of main that performs sensitivity analysis, results of
# which are showin in Extended Data Fig. 1.
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

# data preparation scripts - only these necessary for sensitivity analysis
source("R/00_configure_run.R")
source("R/01_prepare_monthly_data.R")
source("R/02_intercalibrate_periods.R")
source("R/03_set_local_baseline_range.R")
source("R/04_detect_local_deviations.R")
source("R/05_build_ice_areas.R")
source("R/06_get_land_area_with_local_deviations.R")

# data export & plotting script - only this necessary for sensitivity analysis
source("R/output01_land_area_with_local_deviations.R")

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

# local bound options
cfg$local_bound_low <- 0.025 # change this to 0.01 for Extended Data Fig. 2e-f
cfg$local_bound_high <- 0.975 # change this to 0.99 for Extended Data Fig. 2e-f

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
remove(cfg, ensemble_member, variables, gcms, ghms, v, i, j, new_config)

# Prepare main analysis data with prescribed configs ---------------------------
# WILL OVERWRITE ALL FILES IN DATA FOLDER

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

        }
stopImplicitCluster()

# Create outputs from sensitivity analysis (Extended Data Fig. 2) --------------

cfg <- configure_run("configs/dis main isimip2")
output_land_area_with_local_deviations(cfg)

cfg <- configure_run("configs/rootmoist main isimip2")
output_land_area_with_local_deviations(cfg)




