## freshwater-pb

### Code repository for: "Streamflow and soil moisture shift far beyond pre-industrial conditions globally – planetary boundary for freshwater change transgressed"

Miina Porkka, Vili Virkki, Lan Wang-Erlandsson, Dieter Gerten, Tom Gleeson,
Chinchu Mohan, Ingo Fetzer, Fernando Jaramillo, Arie Staal, Sofie te Wierik,
Arne Tobian, Ruud van der Ent, Petra Döll, Martina Flörke, Simon N. Gosling,
Naota Hanasaki, Yusuke Satoh, Hannes Müller Schmied, Niko Wanders,
James S. Famiglietti, Johan Rockström, Matti Kummu

**Published in _Nature Water_**  
**DOI:** <span style="color:red">**ADD LINK WHEN DOI AVAILABLE**</span>.

**Output data available in a Zenodo repository**  
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10531807.svg)](https://doi.org/10.5281/zenodo.10531807)

**Corresponding authors of the article**  
Miina Porkka (miina.porkka@uef.fi)  
Vili Virkki (vili.virkki@aalto.fi)

**Repository author**  
Vili Virkki (vili.virkki@aalto.fi)

**Please cite the _Nature Water_ publication if using code from this repository
or data from the Zenodo repository in another publication.**

### Repository structure

```
freshwater-pb
├── Data
│   ├── anthromes
│   │   └── readme_anthromes.txt
│   └── ardiv
│       ├── hybas2.gpkg
│       ├── hybas2.tif
│       ├── hybas3.gpkg
│       ├── hybas3.tif
│       └── selected_basins.gpkg
├── LICENSE.md
├── R
│   ├── 00_configure_run.R
│   ├── 01_prepare_monthly_data.R
│   ├── 02_intercalibrate_periods.R
│   ├── 03_set_local_baseline_range.R
│   ├── 04_detect_local_deviations.R
│   ├── 05_build_ice_areas.R
│   ├── 06_get_land_area_with_local_deviations.R
│   ├── 07_get_local_deviation_frequency.R
│   ├── isimip_data_download.R
│   ├── main.R
│   ├── output01_land_area_with_local_deviations.R
│   ├── output02_local_deviation_frequency_change.R
│   ├── plot01_regional_persistent_transgressions.R
│   ├── plot02_combined_local_deviation_frequency_increases.R
│   └── sensitivity_analysis.R
├── README.md
├── configs
│   ├── dis main isimip2
│   ├── dis visopts
│   ├── rootmoist main isimip2
│   ├── rootmoist visopts
│   └── template main isimip2
├── logs
│   ├── log files – not listed explicitly here
│   └── R sessionInfo() output
└── txt
    ├── complete_data_folder_tree.txt
    ├── dis_raw_data_filelist.txt
    ├── dis_raw_data_tree.txt
    ├── rootmoist_raw_data_filelist.txt
    └── rootmoist_raw_data_tree.txt

```

### Folders

The repository is divided into the following folders.

#### Data

Upon execution, the R scripts read from `Data/ardiv` and `Data/anthromes` and thus
these folders should not be altered. Otherwise, the scripts will create their
required folder structure and add intermediate outputs within `Data`. Final outputs
and graphics are written into folder `output` that is created in repository root
when an output script is called.

Raw hydrological data is not given with this repository, but information on
acquiring raw data is given below. Output data is available in a separate repository
listed above.

This code repository contains data for two levels of rasterised
[HydroBASINS](https://www.hydrosheds.org/products/hydrobasins),
which are used in the manuscript Figure 4 and Figure 5. Additionally, a folder
for [HYDE 3.2](https://essd.copernicus.org/articles/9/927/2017/essd-9-927-2017.html)
anthromes is initialised. Anthromes are required to run the analysis, and users
should download them from the above source.

#### R

`00_` to `07_`: prepare data from raw ISIMIP NetCDFs and perform analysis following
the methodology outlined in Fig. 1.

`output01_`: plot and write output csv files for the land area
with local deviations (used in Fig. 2, 5, Extended Data Fig. 1, 2, 3, 4)

`output02_`: plot and write output csv files for local deviation frequency changes
(used in Fig. 3, Extended Data Fig. 5, 6, 7).

`plot01_`: plot regional time of persistent transgression (used in Fig. 4).

`plot02_`: combine and plot streamflow and soil moisture local deviation frequency
increases (used in Fig. 5, Extended Data Fig. 7).

`main.R` initialises all functions within data preparation scripts and calls them
in order to reproduce analysis and outputs shown in the manuscript.

* Properties common for both `main.R` and `sensitivity_analysis.R`:

    * assume working directory in the repository root

    * be prepared so that the script can run directly from command line
    
    * silently overwrite intermediate outputs in `Data` (make sure to have backups after each run)

    * call required libraries within; tested versions are specified in R `sessionInfo()` output given in `logs`

`sensitivity_analysis.R` initialises functions and performs the workflow needed
for parts of Extended Data Figure 2 in the manuscript.

`isimip_data_download.R` is a helper script for downloading raw data - more of
this below.

#### configs

This folder contains plaintext files that describe parameters and visual options
used throughout the analysis. The format of the files follows a two-piece convention
(separated by semicolon) in which the first piece denotes the parameter name and
the second piece denotes the parameter value.

The config files describing analysis parameters (tagged with `main`) are read
explicitly by `R/00_configure_run.R` and a list of parameters is then passed to
data preparation (`01_` to `07_`) and output (`output01_` & `output02_`) scripts.
Plotting scripts (`plot01_` & `plot02_`) do not need a parameter list but assume
that output scripts have already been run.

The config files describing visual options (tagged with `visopts`) are read only internally
within `output01_` and `output02_` scripts and not directly from `main.R`.

#### logs

The data preparation and output scripts do some logging to plaintext files that
are initialised by `R/00_configure_run.R`. This folder contains log files from the
analysis run from which the outputs shown in the manuscript and shared in the
separate data repository were prepared.

Additionally, `logs` contains an R `sessionInfo()` output that describes which
R and package versions were used in the analysis.

#### txt

This folder contains lists of raw data files and example file trees:

`complete_data_folder_tree.txt`: tree of the complete `Data` folder after running
`main.R`.

`_raw_data_filelist.txt`: URLs that were used in downloading data from the ISIMIP
repository.

* note: a MATSIRO-HadGEM2 file is missing from `dis_raw_data_filelist.txt` (see below)

`_raw_data_tree.txt`: trees of raw data organisation within a directory specified
by config parapeter `raw_data_root`.

### Raw data acquisition

Raw hydrological data should be acquired from [ISIMIP](https://data.isimip.org/search/)
and placed into a directory with the following structure:

```
CONFIG_RAW_DATA_ROOT
└── CONFIG_VARIABLE
    ├── CONFIG_PERIOD_LABEL_A
    │   └── raw
    │       └── add files here
    └── CONFIG_PERIOD_LABEL_B
        └── raw
            └── add files here
```

For the analysis presented in the manuscript, the raw data was stored as shown in
`txt/dis_raw_data_tree.txt` and `txt/rootmoist_raw_data_tree.txt`. A helper script
for downloading ISIMIP data is available in `R/isimip_data_download.R`.

**NOTE:** At the time of downloading ISIMIP 2b data for this analysis, one NetCDF
file for MATSIRO-HadGEM2 (matsiro_hadgem2-es_ewembi_picontrol_1860soc_co2_dis_global_daily_1661_1670.nc4)
was missing from the ISIMIP repository. This file can be requested from
[ISIMIP MATSIRO representatives](https://www.isimip.org/impactmodels/details/91/).

### Changelog and versioning

Should the repository be updated with new versions, main changes will be briefly
summarised here. Commits describing new versions are marked respectively with
Git tags.

#### v1.0.0.

Version used in producing the results shown in the
<span style="color:red">**ADD LINK WHEN DOI AVAILABLE**</span>.
[_Nature Water_ manuscript](https://www.nature.com/natwater).

### License

Attribution 4.0 International (CC BY 4.0)  
https://creativecommons.org/licenses/by/4.0/

