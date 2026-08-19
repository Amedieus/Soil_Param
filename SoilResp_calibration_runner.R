library(PEcAn.all)
library(PEcAn.SIPNET)
library(PEcAn.uncertainty)
library(PEcAn.settings)
library(future)
library(furrr)
library(data.table)
library(neonSoilFlux)
source("/projectnb/dietzelab/guYANG/soilparam/SoilResp_calibration_function.R")
source("/projectnb/dietzelab/guYANG/soilparam/SoilT_tau_functions.R")
source("/projectnb/dietzelab/guYANG/soilparam/SoilT_tau_permafrost_functions.R")
setwd("/projectnb/dietzelab/guYANG/soilparam")

###Data prep
## index
lookup <- fread("/projectnb/dietzelab/guYANG/soilparam/lookup.csv")
## Get NEON soilresp observations for multi site
multi_site <- get_soil_neon_data_multi_site(
  lookup = lookup,
  start_date = "2017-01-01",
  end_date = "2024-12-31",
  output_dir = "/projectnb/dietzelab/guYANG/soilparam/NEON_calibration_data"
)
save(multi_site, file = "/projectnb/dietzelab/guYANG/soilparam/soilphysic_NEON.RData")

## Get NEON sites' MET drivers
era5_all <- get_soil_era5_data_multi_site(
  lookup = lookup,
  start_date = "2017-01-01",
  end_date = "2024-12-31",
  output_dir = "/projectnb/dietzelab/guYANG/soilparam/calibration_data"
)
save(era5_all, file = "/projectnb/dietzelab/guYANG/soilparam/ear5_NEON.RData")

### Optimized tau-MLE for non-permafrost
non_permafrost_tau <- run_all_indices_depths_parallel(
  lookup = non_permafrost,
  workers = 16L,
  obs_start_year = 2017L,
  obs_end_year = 2024L
)

### Optimized tau-MLE for permafrost (Alaska)
permafrost_tau_results <- fit_permafrost_all_sites(
  lookup = permafrost,
  start_year = 2017L,
  end_year = 2024L,
  multi_site = multi_site,
  era5_all = era5_all,
  vertical_positions = "502",
  workers = 5L,
  # Tau-only MLE
  fit_tau_only = TRUE,
  fixed_a_C = 0,
  fixed_n_warm = 1,
  fixed_n_cold = 0.3,
  tau_bounds = c(
    0.125,
    180
  ),
  # Data / LOYO
  warmup_days = 180L,
  min_observations = 100L,
  min_days_per_year = 120L,
  min_days_per_season = 10L,
  output_dir =
    "/projectnb/dietzelab/guYANG/soilparam/permafrost_tau_only"
)

### SM validation
soilmoisture_compare <- prepare_soilmoisture_validation_table(
  lookup = lookup,
  multi_site = multi_site,
  era5_all = era5_all,
  start_date = "2017-01-01",
  end_date = "2024-12-31",
  sipnet_out_dir =
    "/projectnb/dietzelab/guYANG/pecan/updated_clim/out",
  output_file =
    "/projectnb/dietzelab/guYANG/soilparam/SoilMoisture_validation/NEON_SIPNET_SoilMoisture_3hour.csv"
)