# ------------------------------------------------------------------------------
# 
# isimip_data_download
# 
# Download data from ISIMIP data portal using a list of URLs, with download
# directory specified manually and assumed folder structure prepared.
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

library(curl)
library(doParallel)
library(stringr)
library(dplyr)

# replace in out_loc:
# CONFIG_RAW_DATA_ROOT = root folder for downloading data, same as in config file
# CONFIG_VARIABLE = which variable is downloaded, here dis or rootmoist

# required folder structure
# CONFIG_RAW_DATA_ROOT/
# └── CONFIG_VARIABLE
#     ├── historical
#     │   └── raw
#     └── preind
#         └── raw

# filelists (dis & rootmoist) are in public repository
files <- readLines("dis_raw_data_filelist.txt")

registerDoParallel(cores = 1)
foreach(i = 1:length(files),
        .errorhandling = "remove",
        .verbose = TRUE) %dopar% {
          out <- files[i] %>%
            str_split("/") %>%
            unlist() %>%
            tail(1)
          period <- ifelse(grepl("historical", out), "/historical/raw/", "/preind/raw/")
          out_loc <- paste0("/CONFIG_RAW_DATA_ROOT/CONFIG_VARIABLE", period, out)
          # out_loc <- paste0("/Volumes/T7/rootmoist", period, out) # e.g. like this
          h <- new_handle()
          handle_setopt(h, ssl_verifyhost = 0, ssl_verifypeer = 0)
          if (!file.exists(out_loc)){
            curl_download(files[i], out_loc, handle = h)
          }
        }
stopImplicitCluster()
