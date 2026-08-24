suppressPackageStartupMessages({
  library(data.table)
})
source(
  "/projectnb/dietzelab/guYANG/Soil_Param/generate_new_clim_functions.R"
)
INPUT_ROOT <-
  "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024"
OUTPUT_ROOT <-
  "/projectnb/dietzelab/guYANG/pecan/modified_met/new_calibrated_clim"

NON_PERMAFROST_OUTPUT <-
  file.path(OUTPUT_ROOT)

PERMAFROST_OUTPUT <-
  file.path(OUTPUT_ROOT)

MEMBERS <- 1:10
START_DATE <- "2012-01-01"
END_DATE <- "2024-12-31"
N_CORES <- 18L

dir.create(
  NON_PERMAFROST_OUTPUT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  PERMAFROST_OUTPUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. 准备两个 lookup
# ============================================================

non_permafrost_lookup <-
  as.data.table(
    copy(non_permafrost)
  )

permafrost_lookup <-
  as.data.table(
    copy(permafrost)
  )

required_lookup_columns <- c(
  "index",
  "final_pft"
)

if (
  !all(
    required_lookup_columns %in%
    names(non_permafrost_lookup)
  )
) {
  stop(
    "non_permafrost must contain index and final_pft."
  )
}

if (
  !all(
    required_lookup_columns %in%
    names(permafrost_lookup)
  )
) {
  stop(
    "permafrost must contain index and final_pft."
  )
}

if (anyDuplicated(non_permafrost_lookup$index)) {
  stop(
    "non_permafrost contains duplicated index values."
  )
}

if (anyDuplicated(permafrost_lookup$index)) {
  stop(
    "permafrost contains duplicated index values."
  )
}

overlapping_indices <- intersect(
  non_permafrost_lookup$index,
  permafrost_lookup$index
)

if (length(overlapping_indices) > 0L) {
  stop(
    "The same indices occur in both lookups: ",
    paste(
      head(overlapping_indices, 20L),
      collapse = ", "
    )
  )
}


# ============================================================
# 3. non-permafrost PFT tau 参数
# ============================================================

# pft_tau_test$pft_tau 应该一行对应一个 final_pft，并包含：
# final_pft
# pft_tau_days

non_permafrost_pft_parameters <-
  as.data.table(
    copy(
      pft_tau_test$pft_tau
    )
  )

required_non_permafrost_columns <- c(
  "final_pft",
  "pft_tau_days"
)

if (
  !all(
    required_non_permafrost_columns %in%
    names(non_permafrost_pft_parameters)
  )
) {
  stop(
    "pft_tau_test$pft_tau must contain final_pft and pft_tau_days."
  )
}

if (
  anyDuplicated(
    non_permafrost_pft_parameters$final_pft
  )
) {
  stop(
    "The non-permafrost PFT parameter table has duplicated final_pft."
  )
}


# ============================================================
# 4. 将 permafrost site 参数聚合到 PFT
# ============================================================

# permafrost_tau_results$summary 应包含：
# index, a_C, n_warm, n_cold, tau_days

permafrost_site_parameters <-
  as.data.table(
    copy(
      permafrost_tau_results$summary
    )
  )

# 如果结果中包含多个深度，只使用 VER501
if (
  "verticalPosition" %in%
  names(permafrost_site_parameters)
) {
  permafrost_site_parameters <-
    permafrost_site_parameters[
      as.character(verticalPosition) == "501"
    ]
}

if (
  "vertical_position" %in%
  names(permafrost_site_parameters)
) {
  permafrost_site_parameters <-
    permafrost_site_parameters[
      as.character(vertical_position) == "501"
    ]
}

required_permafrost_columns <- c(
  "index",
  "a_C",
  "n_warm",
  "n_cold",
  "tau_days"
)

if (
  !all(
    required_permafrost_columns %in%
    names(permafrost_site_parameters)
  )
) {
  stop(
    paste0(
      "permafrost_tau_results$summary must contain: ",
      paste(
        required_permafrost_columns,
        collapse = ", "
      )
    )
  )
}

permafrost_index_to_pft <-
  unique(
    permafrost_lookup[
      ,
      .(
        index,
        final_pft
      )
    ]
  )

permafrost_site_parameters <-
  merge(
    permafrost_site_parameters[
      ,
      .(
        index,
        a_C,
        n_warm,
        n_cold,
        tau_days
      )
    ],
    permafrost_index_to_pft,
    by = "index",
    all = FALSE
  )

permafrost_site_parameters <-
  permafrost_site_parameters[
    is.finite(a_C) &
      is.finite(n_warm) &
      is.finite(n_cold) &
      is.finite(tau_days) &
      n_warm >= 0 &
      n_cold >= 0 &
      tau_days > 0
  ]

if (nrow(permafrost_site_parameters) == 0L) {
  stop(
    "No valid permafrost site parameter rows remain."
  )
}

# 与普通 PFT tau workflow 一致，这里使用 PFT 内 site MLE 的均值
permafrost_pft_parameters <-
  permafrost_site_parameters[
    ,
    .(
      a_C =
        mean(
          a_C,
          na.rm = TRUE
        ),
      
      n_warm =
        mean(
          n_warm,
          na.rm = TRUE
        ),
      
      n_cold =
        mean(
          n_cold,
          na.rm = TRUE
        ),
      
      pft_tau_days =
        mean(
          tau_days,
          na.rm = TRUE
        ),
      
      n_calibration_sites =
        uniqueN(index)
    ),
    by =
      final_pft
  ]

setorder(
  permafrost_pft_parameters,
  final_pft
)

fwrite(
  permafrost_pft_parameters,
  file.path(
    PERMAFROST_OUTPUT,
    "permafrost_pft_process_parameters.csv"
  )
)

print(
  permafrost_pft_parameters
)


# ============================================================
# 5. 生成 non-permafrost clim
# ============================================================

non_permafrost_manifest <-
  generate_tau_clims_for_all_pfts(
    lookup =
      non_permafrost_lookup,
    
    newpft =
      non_permafrost_lookup,
    
    pft_tau_test =
      non_permafrost_pft_parameters,
    
    input_root =
      INPUT_ROOT,
    
    output_root =
      NON_PERMAFROST_OUTPUT,
    
    members =
      MEMBERS,
    
    start_date =
      START_DATE,
    
    end_date =
      END_DATE,
    
    n_cores =
      N_CORES,
    
    overwrite =
      FALSE,
    
    verbose =
      TRUE,
    
    missing_tau_action =
      "error",
    
    soil_temperature_model =
      "non_permafrost",
    
    clim_format_version =
      "v2",
    
    clamp_soil_vpd =
      TRUE,
    
    stop_on_error =
      FALSE
  )


# ============================================================
# 6. 生成 permafrost clim
# ============================================================

permafrost_manifest <-
  generate_tau_clims_for_all_pfts(
    lookup =
      permafrost_lookup,
    
    newpft =
      permafrost_lookup,
    
    pft_tau_test =
      permafrost_pft_parameters,
    
    input_root =
      INPUT_ROOT,
    
    output_root =
      PERMAFROST_OUTPUT,
    
    members =
      MEMBERS,
    
    start_date =
      START_DATE,
    
    end_date =
      END_DATE,
    
    n_cores =
      N_CORES,
    
    overwrite =
      FALSE,
    
    verbose =
      TRUE,
    
    missing_tau_action =
      "error",
    
    soil_temperature_model =
      "permafrost",
    
    clim_format_version =
      "v2",
    
    clamp_soil_vpd =
      TRUE,
    
    stop_on_error =
      FALSE
  )


# ============================================================
# 7. 合并并检查 manifest
# ============================================================

non_permafrost_manifest <-
  as.data.table(
    non_permafrost_manifest
  )

permafrost_manifest <-
  as.data.table(
    permafrost_manifest
  )

combined_manifest <-
  rbindlist(
    list(
      non_permafrost_manifest,
      permafrost_manifest
    ),
    use.names = TRUE,
    fill = TRUE
  )

fwrite(
  combined_manifest,
  file.path(
    OUTPUT_ROOT,
    "combined_new_clim_manifest.csv"
  )
)

cat(
  "\n==============================\n",
  "NON-PERMAFROST STATUS\n",
  "==============================\n"
)

print(
  non_permafrost_manifest[
    ,
    .N,
    by = status
  ]
)

cat(
  "\n==============================\n",
  "PERMAFROST STATUS\n",
  "==============================\n"
)

print(
  permafrost_manifest[
    ,
    .N,
    by = status
  ]
)

failed_jobs <-
  combined_manifest[
    status %in%
      c(
        "MISSING_INPUT",
        "ERROR",
        "PFT_ERROR"
      )
  ]

cat(
  "\nFailed jobs:",
  nrow(failed_jobs),
  "\n"
)

if (nrow(failed_jobs) > 0L) {
  fwrite(
    failed_jobs,
    file.path(
      OUTPUT_ROOT,
      "failed_new_clim_jobs.csv"
    )
  )
}

cat(
  "\nNew climate files are under:\n",
  normalizePath(
    OUTPUT_ROOT,
    mustWork = FALSE
  ),
  "\n"
)