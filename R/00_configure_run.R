# ------------------------------------------------------------------------------
# 
# 00_configure_run
# 
# Read a configuration file describing a run setup of GHMs and GCMs together
# with parameters on how the analysis should be run.
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

configure_run <- function(config_path) {

  # read config file; return a list of options to be used in analysis
  cfg_data <- readLines(config_path) %>%
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

  cfg$tempdir <- tempdir()
  cfg$is_corrected <- FALSE
  if (!dir.exists("logs")) { dir.create("logs") }
  cfg$log_file <- paste0("logs/log ", Sys.time(), " ", cfg$variable, " ",
                         paste0(paste(cfg$impactmodel, collapse = "_"), " ",
                                paste(cfg$forcing, collapse = "_"), ".txt"))

  write(format(Sys.time(), "%a %b %d %X %Y"), cfg$log_file, append = TRUE)
  write("opening log...", cfg$log_file, append = TRUE)

  return (cfg)

}
