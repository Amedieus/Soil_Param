library(PEcAn.all)
library(PEcAn.SIPNET)
library(PEcAn.uncertainty)
library(PEcAn.settings)
library(future)
library(furrr)
library(data.table)
library(neonSoilFlux)

project_root <- "/projectnb/dietzelab/guYANG/soilparam"

source(file.path(project_root, "SoilResp_calibration_function.R"))
source(file.path(project_root, "neon_tsoil_mle_workflow.R"))

setwd(project_root)

# =============================================================================
# Data preparation
# =============================================================================

lookup <- data.table::fread(
  file.path(project_root, "lookup.csv")
)

multi_site <- get_soil_neon_data_multi_site(
  lookup = lookup,
  start_date = "2017-01-01",
  end_date = "2024-12-31",
  output_dir = file.path(project_root, "NEON_calibration_data")
)

save(
  multi_site,
  file = file.path(project_root, "soilphysic_NEON.RData")
)

era5_all <- get_soil_era5_data_multi_site(
  lookup = lookup,
  start_date = "2017-01-01",
  end_date = "2024-12-31",
  output_dir = file.path(project_root, "calibration_data")
)

save(
  era5_all,
  file = file.path(project_root, "era5_NEON.RData")
)

# =============================================================================
# Unified site-specific SoilT MLE + LOYO
#
# The function performs the latitude classification internally:
#   site_lat > 63  -> permafrost
#   site_lat <= 63 -> non-permafrost
# =============================================================================

neon_mle_results <- neon_tsoil_mle(
  lookup = lookup,
  start_year = 2017L,
  end_year = 2024L,
  multi_site = multi_site,
  era5_all = era5_all,
  vertical_positions = "502",
  workers = 16L,
  permafrost_latitude_threshold = 63,
  permafrost_fit_tau_only = TRUE,
  permafrost_fixed_a_C = 0,
  permafrost_fixed_n_warm = 1,
  permafrost_fixed_n_cold = 0.3,
  tau_bounds = c(0.125, 180),
  warmup_days = 180L,
  min_observations = 100L,
  min_days_per_year = 120L,
  min_days_per_season = 10L,
  output_dir = file.path(project_root, "neon_tsoil_mle")
)

# =============================================================================
# Soil-moisture validation
# =============================================================================

soilmoisture_compare <- prepare_soilmoisture_validation_table(
  lookup = lookup,
  multi_site = multi_site,
  era5_all = era5_all,
  start_date = "2017-01-01",
  end_date = "2024-12-31",
  sipnet_out_dir =
    "/projectnb/dietzelab/guYANG/pecan/updated_clim/out",
  output_file = file.path(
    project_root,
    "SoilMoisture_validation",
    "NEON_SIPNET_SoilMoisture_3hour.csv"
  )
)

