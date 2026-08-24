# =============================================================================
# Generate SIPNET climate files using NEON SoilT first and site tau for gaps
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# =============================================================================
# 1. User configuration
# =============================================================================

target_vertical_position <- "502"
start_date <- "2012-01-01"
end_date <- "2024-12-31"

soil_param_code_dir <-
  "/projectnb/dietzelab/guYANG/Soil_Param"

lookup_file <-
  "/projectnb/dietzelab/guYANG/soilparam/lookup.csv"

multi_site_file <-
  "/projectnb/dietzelab/guYANG/soilparam/soilphysic_NEON.RData"

non_permafrost_parameter_file <- file.path(
  "/projectnb/dietzelab/guYANG/soilparam/tau_by_index",
  "all_indices_all_depths_summary.csv"
)

permafrost_parameter_file <- file.path(
  "/projectnb/dietzelab/guYANG/soilparam/permafrost_tau_only",
  "permafrost_tsoil_tau_only_MLE_summary.csv"
)

input_clim_root <-
  "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024"

output_clim_root <- file.path(
  soil_param_code_dir,
  paste0(
    "neon_tau_gapfilled_clim_VER",
    target_vertical_position
  )
)


# =============================================================================
# 2. Load functions and inputs
# =============================================================================

source(
  file.path(
    soil_param_code_dir,
    "generate_new_clim_functions.R"
  )
)

required_files <- c(
  lookup_file,
  multi_site_file,
  non_permafrost_parameter_file,
  permafrost_parameter_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

lookup <- data.table::fread(
  lookup_file
)

load(
  multi_site_file
)

if (!exists("multi_site", inherits = FALSE)) {
  stop(
    "The multi-site RData file did not create an object named `multi_site`.",
    call. = FALSE
  )
}

non_permafrost_parameters <- data.table::fread(
  non_permafrost_parameter_file
)

permafrost_parameters <- data.table::fread(
  permafrost_parameter_file
)


# =============================================================================
# 3. Restrict lookup to sites with calibrated tau parameters
# =============================================================================

calibrated_indices <- sort(
  unique(
    c(
      as.integer(non_permafrost_parameters$index),
      as.integer(permafrost_parameters$index)
    )
  )
)

calibration_lookup <- lookup[
  as.integer(index) %in% calibrated_indices
]

if (nrow(calibration_lookup) == 0L) {
  stop(
    "No lookup rows match the calibrated tau parameter tables.",
    call. = FALSE
  )
}


# =============================================================================
# 4. Generate updated SIPNET climate files
# =============================================================================
# 明确只保留502
non_permafrost_parameters <- non_permafrost_parameters[
  as.character(vertical_position) ==
    as.character(target_vertical_position)
]

# 应该每个index只剩一行
non_permafrost_parameters[
  ,
  .N,
  by = index
][
  N != 1L
]


neon_tau_clim_manifest <- generate_neon_tau_gapfilled_clims(
  lookup = calibration_lookup,
  multi_site = multi_site,
  
  non_permafrost_parameters =
    non_permafrost_parameters,
  
  permafrost_parameters =
    permafrost_parameters,
  
  vertical_position =
    target_vertical_position,
  
  # 真正控制新文件的时间范围
  start_date = "2012-01-01",
  end_date = "2024-12-31",
  
  input_root = input_clim_root,
  output_root = output_clim_root,
  
  # 生成10个集合成员
  members = 1:10,
  
  # 只用于定位原始clim文件
  input_clim_start_date = "2012-01-01",
  input_clim_end_date = "2024-12-31",
  
  permafrost_latitude_threshold = 63,
  
  permafrost_a_C = 0,
  permafrost_n_warm = 1,
  permafrost_n_cold = 0.3,
  
  workers = 1L,
  warmup_days = 180L,
  
  observation_match_tolerance_seconds = 60,
  
  overwrite = FALSE,
  write_diagnostics = TRUE,
  stop_on_error = FALSE
)