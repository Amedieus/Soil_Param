prepare_soilmoisture_validation_table <-
function(
    lookup,
    multi_site,
    era5_all,
    start_date,
    end_date,
    sipnet_out_dir =
      "/projectnb/dietzelab/guYANG/pecan/updated_clim/out",
    output_file =
      "/projectnb/dietzelab/guYANG/soilparam/SoilMoisture_validation/NEON_SIPNET_SoilMoisture_3hour.csv",
    keep_only_paired = TRUE,
    include_era5_drivers = TRUE
) {
  
  # ==========================================================================
  # 1. Validate prepared input objects
  # ==========================================================================
  
  neon_sm <- getElement(
    multi_site,
    "soilMoisture"
  )
  
  
  era5_mean <- getElement(
    era5_all,
    "era5_mean"
  )
  
  
  if (
    is.null(neon_sm) ||
    nrow(neon_sm) == 0L
  ) {
    
    stop(
      "`multi_site` does not contain valid `soilMoisture` data.",
      call. = FALSE
    )
  }
  
  
  if (
    is.null(era5_mean) ||
    nrow(era5_mean) == 0L
  ) {
    
    stop(
      "`era5_all` does not contain valid `era5_mean` data.",
      call. = FALSE
    )
  }
  
  
  neon_sm <- data.table::as.data.table(
    data.table::copy(
      neon_sm
    )
  )
  
  
  era5_mean <- data.table::as.data.table(
    data.table::copy(
      era5_mean
    )
  )
  
  
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
  
  
  # ==========================================================================
  # 2. Validate required columns
  # ==========================================================================
  
  required_neon_columns <- c(
    "time",
    "index",
    "verticalPosition",
    "VSWC"
  )
  
  
  missing_neon_columns <- setdiff(
    required_neon_columns,
    names(
      neon_sm
    )
  )
  
  
  if (length(missing_neon_columns) > 0L) {
    
    stop(
      "`multi_site$soilMoisture` is missing: ",
      paste(
        missing_neon_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  required_era5_columns <- c(
    "time",
    "index"
  )
  
  
  missing_era5_columns <- setdiff(
    required_era5_columns,
    names(
      era5_mean
    )
  )
  
  
  if (length(missing_era5_columns) > 0L) {
    
    stop(
      "`era5_all$era5_mean` is missing: ",
      paste(
        missing_era5_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  
  if (
    !"index" %in%
    names(
      lookup_dt
    )
  ) {
    
    stop(
      "`lookup` must contain `index`.",
      call. = FALSE
    )
  }
  
  
  # ==========================================================================
  # 3. Requested time range
  # ==========================================================================
  
  start_time <- as.POSIXct(
    paste0(
      start_date,
      " 00:00:00"
    ),
    tz = "UTC"
  )
  
  
  end_time <- as.POSIXct(
    paste0(
      end_date,
      " 00:00:00"
    ),
    tz = "UTC"
  ) +
    86400
  
  
  if (
    is.na(start_time) ||
    is.na(end_time)
  ) {
    
    stop(
      "Dates must use YYYY-MM-DD.",
      call. = FALSE
    )
  }
  
  
  if (end_time <= start_time) {
    
    stop(
      "`end_date` must be >= `start_date`.",
      call. = FALSE
    )
  }
  
  
  # ==========================================================================
  # 4. Determine requested model indices
  # ==========================================================================
  
  lookup_indices <- unique(
    as.integer(
      lookup_dt[
        !is.na(index),
        index
      ]
    )
  )
  
  
  lookup_indices <- lookup_indices[
    is.finite(
      lookup_indices
    )
  ]
  
  
  neon_indices <- unique(
    as.integer(
      neon_sm$index
    )
  )
  
  
  era5_indices <- unique(
    as.integer(
      era5_mean$index
    )
  )
  
  
  target_indices <- sort(
    Reduce(
      intersect,
      list(
        lookup_indices,
        neon_indices,
        era5_indices
      )
    )
  )
  
  
  if (length(target_indices) == 0L) {
    
    stop(
      paste0(
        "No common model indices were found among lookup, ",
        "multi_site, and era5_all."
      ),
      call. = FALSE
    )
  }
  
  
  # ==========================================================================
  # 5. Prepare NEON soil moisture
  # ==========================================================================
  
  neon_sm[
    ,
    time :=
      as.POSIXct(
        time,
        tz = "UTC"
      )
  ]
  
  
  neon_sm[
    ,
    index :=
      as.integer(
        index
      )
  ]
  
  
  neon_sm[
    ,
    verticalPosition :=
      as.character(
        verticalPosition
      )
  ]
  
  
  neon_sm <- neon_sm[
    index %in%
      target_indices &
      time >=
      start_time &
      time <
      end_time
  ]
  
  
  data.table::setnames(
    neon_sm,
    "VSWC",
    "VSWC_NEON"
  )
  
  
  if (
    "VSWC_temporal_sd" %in%
    names(
      neon_sm
    )
  ) {
    
    data.table::setnames(
      neon_sm,
      "VSWC_temporal_sd",
      "VSWC_NEON_temporal_sd"
    )
  }
  
  
  if (
    "VSWC_spatial_sd" %in%
    names(
      neon_sm
    )
  ) {
    
    data.table::setnames(
      neon_sm,
      "VSWC_spatial_sd",
      "VSWC_NEON_spatial_sd"
    )
  }
  
  
  # ==========================================================================
  # 6. Prepare ERA5 exact-time meteorological table
  # ==========================================================================
  
  era5_mean[
    ,
    time :=
      as.POSIXct(
        time,
        tz = "UTC"
      )
  ]
  
  
  era5_mean[
    ,
    index :=
      as.integer(
        index
      )
  ]
  
  
  era5_mean <- era5_mean[
    index %in%
      target_indices &
      time >=
      start_time &
      time <
      end_time
  ]
  
  
  # --------------------------------------------------------------------------
  # Keep only one ERA5 record per index and exact timestamp
  # --------------------------------------------------------------------------
  
  era5_mean <- unique(
    era5_mean,
    by = c(
      "index",
      "time"
    )
  )
  
  
  # ==========================================================================
  # 7. Helper functions for SIPNET extraction
  # ==========================================================================
  
  mean_finite <- function(
    x
  ) {
    
    x <- suppressWarnings(
      as.numeric(
        x
      )
    )
    
    
    x <- x[
      is.finite(
        x
      )
    ]
    
    
    if (length(x) == 0L) {
      
      return(
        NA_real_
      )
    }
    
    
    return(
      mean(
        x
      )
    )
  }
  
  
  sd_finite <- function(
    x
  ) {
    
    x <- suppressWarnings(
      as.numeric(
        x
      )
    )
    
    
    x <- x[
      is.finite(
        x
      )
    ]
    
    
    if (length(x) < 2L) {
      
      return(
        NA_real_
      )
    }
    
    
    return(
      stats::sd(
        x
      )
    )
  }
  
  
  # ==========================================================================
  # 8. Find SIPNET ensemble directories once
  # ==========================================================================
  
  all_ensemble_dirs <- list.dirs(
    sipnet_out_dir,
    recursive = FALSE,
    full.names = TRUE
  )
  
  
  if (length(all_ensemble_dirs) == 0L) {
    
    stop(
      "No SIPNET ensemble directories were found in: ",
      sipnet_out_dir,
      call. = FALSE
    )
  }
  
  
  # ==========================================================================
  # 9. Extract SIPNET soil moisture for one index
  # ==========================================================================
  
  read_sipnet_index <- function(
    index_i
  ) {
    
    directory_pattern <- paste0(
      "^ENS-[0-9]+-",
      index_i,
      "$"
    )
    
    
    ensemble_dirs <- all_ensemble_dirs[
      grepl(
        directory_pattern,
        basename(
          all_ensemble_dirs
        )
      )
    ]
    
    
    if (length(ensemble_dirs) == 0L) {
      
      return(
        NULL
      )
    }
    
    
    ensemble_data <- vector(
      "list",
      length(
        ensemble_dirs
      )
    )
    
    
    for (
      ensemble_i in
      seq_along(
        ensemble_dirs
      )
    ) {
      
      ensemble_dir <- ensemble_dirs[
        ensemble_i
      ]
      
      
      sipnet_file <- file.path(
        ensemble_dir,
        "sipnet.out"
      )
      
      
      if (!file.exists(sipnet_file)) {
        
        next
      }
      
      
      # -----------------------------------------------------------------------
      # Read header to identify SIPNET columns
      # -----------------------------------------------------------------------
      
      header <- tryCatch(
        names(
          data.table::fread(
            sipnet_file,
            nrows = 0L,
            showProgress = FALSE
          )
        ),
        error = function(e) {
          character()
        }
      )
      
      
      if (length(header) == 0L) {
        
        next
      }
      
      
      header_lower <- tolower(
        header
      )
      
      
      year_column <- header[
        header_lower ==
          "year"
      ]
      
      
      day_column <- header[
        header_lower %in%
          c(
            "day",
            "doy"
          )
      ]
      
      
      time_column <- header[
        header_lower ==
          "time"
      ]
      
      
      soilwater_column <- header[
        header_lower ==
          "soilwater"
      ]
      
      
      wetness_column <- header[
        header_lower ==
          "soilwetnessfrac"
      ]
      
      
      if (
        length(year_column) != 1L ||
        length(day_column) != 1L ||
        length(time_column) != 1L ||
        length(soilwater_column) != 1L ||
        length(wetness_column) != 1L
      ) {
        
        next
      }
      
      
      # -----------------------------------------------------------------------
      # Read only required variables
      # -----------------------------------------------------------------------
      
      x <- tryCatch(
        data.table::fread(
          sipnet_file,
          select = c(
            year_column,
            day_column,
            time_column,
            soilwater_column,
            wetness_column
          ),
          showProgress = FALSE
        ),
        error = function(e) {
          NULL
        }
      )
      
      
      if (
        is.null(x) ||
        nrow(x) == 0L
      ) {
        
        next
      }
      
      
      data.table::setnames(
        x,
        old = c(
          year_column,
          day_column,
          time_column,
          soilwater_column,
          wetness_column
        ),
        new = c(
          "year_raw",
          "doy_raw",
          "hour_raw",
          "soilWater",
          "soilWetnessFrac"
        )
      )
      
      
      x[
        ,
        `:=`(
          year_raw =
            as.integer(
              year_raw
            ),
          
          doy_raw =
            as.numeric(
              doy_raw
            ),
          
          hour_raw =
            as.numeric(
              hour_raw
            ),
          
          soilWater =
            as.numeric(
              soilWater
            ),
          
          soilWetnessFrac =
            as.numeric(
              soilWetnessFrac
            )
        )
      ]
      
      
      # -----------------------------------------------------------------------
      # Reconstruct exact 3-hour SIPNET time
      #
      # SIPNET floating-point time is used only for the first timestamp.
      # Every later timestamp follows exact row order at 3-hour intervals.
      # -----------------------------------------------------------------------
      
      first_hour_exact <- round(
        x$hour_raw[1L] /
          3
      ) *
        3
      
      
      first_datetime <- as.POSIXct(
        sprintf(
          "%04d-01-01 00:00:00",
          x$year_raw[1L]
        ),
        tz = "UTC"
      ) +
        (
          x$doy_raw[1L] -
            1
        ) *
        86400 +
        first_hour_exact *
        3600
      
      
      x[
        ,
        time :=
          first_datetime +
          (
            seq_len(
              .N
            ) -
              1L
          ) *
          10800
      ]
      
      
      x <- x[
        time >=
          start_time &
          time <
          end_time
      ]
      
      
      if (nrow(x) == 0L) {
        
        next
      }
      
      
      ensemble_number <- as.integer(
        sub(
          "^ENS-([0-9]+)-.*$",
          "\\1",
          basename(
            ensemble_dir
          )
        )
      )
      
      
      x[
        ,
        ensemble :=
          ensemble_number
      ]
      
      
      ensemble_data[
        ensemble_i
      ] <- list(
        x[
          ,
          .(
            time,
            ensemble,
            soilWater,
            soilWetnessFrac
          )
        ]
      )
    }
    
    
    valid_ensemble_data <- Filter(
      Negate(
        is.null
      ),
      ensemble_data
    )
    
    
    if (length(valid_ensemble_data) == 0L) {
      
      return(
        NULL
      )
    }
    
    
    all_members <- data.table::rbindlist(
      valid_ensemble_data,
      use.names = TRUE,
      fill = TRUE
    )
    
    
    # -------------------------------------------------------------------------
    # Ensemble mean at every exact 3-hour timestep
    # -------------------------------------------------------------------------
    
    sipnet_mean <- all_members[
      ,
      .(
        SoilWater_SIPNET =
          mean_finite(
            soilWater
          ),
        
        SoilWater_SIPNET_sd =
          sd_finite(
            soilWater
          ),
        
        SoilWetnessFrac_SIPNET =
          mean_finite(
            soilWetnessFrac
          ),
        
        SoilWetnessFrac_SIPNET_sd =
          sd_finite(
            soilWetnessFrac
          ),
        
        n_SIPNET_ensemble =
          data.table::uniqueN(
            ensemble
          )
      ),
      by =
        time
    ]
    
    
    sipnet_mean[
      ,
      index :=
        as.integer(
          index_i
        )
    ]
    
    
    data.table::setorder(
      sipnet_mean,
      time
    )
    
    
    return(
      sipnet_mean
    )
  }
  
  
  # ==========================================================================
  # 10. Extract SIPNET soil moisture for all requested indices
  # ==========================================================================
  
  sipnet_list <- vector(
    "list",
    length(
      target_indices
    )
  )
  
  
  status <- data.table::data.table(
    index =
      target_indices,
    
    sipnet_success =
      FALSE,
    
    n_sipnet_time =
      0L
  )
  
  
  for (
    site_i in
    seq_along(
      target_indices
    )
  ) {
    
    index_i <- target_indices[
      site_i
    ]
    
    
    message(
      "[",
      site_i,
      "/",
      length(
        target_indices
      ),
      "] SIPNET soil moisture | index = ",
      index_i
    )
    
    
    sipnet_i <- read_sipnet_index(
      index_i
    )
    
    
    if (
      is.null(sipnet_i) ||
      nrow(sipnet_i) == 0L
    ) {
      
      next
    }
    
    
    sipnet_list[
      site_i
    ] <- list(
      sipnet_i
    )
    
    
    status[
      site_i,
      `:=`(
        sipnet_success =
          TRUE,
        
        n_sipnet_time =
          nrow(
            sipnet_i
          )
      )
    ]
  }
  
  
  valid_sipnet <- Filter(
    Negate(
      is.null
    ),
    sipnet_list
  )
  
  
  if (length(valid_sipnet) == 0L) {
    
    stop(
      "No valid SIPNET soil-moisture output was found.",
      call. = FALSE
    )
  }
  
  
  sipnet_sm <- data.table::rbindlist(
    valid_sipnet,
    use.names = TRUE,
    fill = TRUE
  )
  
  
  # ==========================================================================
  # 11. Match SIPNET to NEON observations
  #
  # NEON remains the base table because it contains multiple soil depths.
  # ==========================================================================
  
  comparison <- merge(
    neon_sm,
    sipnet_sm,
    by = c(
      "index",
      "time"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  
  # ==========================================================================
  # 12. Add selected ERA5 drivers
  # ==========================================================================
  
  if (
    isTRUE(
      include_era5_drivers
    )
  ) {
    
    era5_driver_columns <- intersect(
      c(
        "index",
        "time",
        "AirT_C",
        "Precip_mm_3h",
        "VPD_Pa",
        "RH",
        "WindSpeed",
        "PAR_umol_m2_s",
        "surface_downwelling_shortwave_flux_in_air",
        "surface_downwelling_longwave_flux_in_air"
      ),
      names(
        era5_mean
      )
    )
    
    
    era5_for_join <- era5_mean[
      ,
      ..era5_driver_columns
    ]
    
    
    comparison <- merge(
      comparison,
      era5_for_join,
      by = c(
        "index",
        "time"
      ),
      all.x = TRUE,
      sort = FALSE
    )
  }
  
  
  # ==========================================================================
  # 13. Add paired-data flag
  # ==========================================================================
  
  comparison[
    ,
    paired :=
      is.finite(
        VSWC_NEON
      ) &
      is.finite(
        SoilWetnessFrac_SIPNET
      )
  ]
  
  
  comparison[
    ,
    n_paired_site_depth :=
      sum(
        paired,
        na.rm = TRUE
      ),
    by = .(
      index,
      verticalPosition
    )
  ]
  
  
  # ==========================================================================
  # 14. Standardize moisture dynamics within each site and depth
  #
  # Raw NEON VSWC and SIPNET soilWetnessFrac are not identical physical
  # quantities, so standardized values are useful for dynamics validation.
  # ==========================================================================
  
  comparison[
    ,
    `:=`(
      VSWC_NEON_z =
        NA_real_,
      
      SoilWetnessFrac_SIPNET_z =
        NA_real_
    )
  ]
  
  
  comparison[
    paired == TRUE,
    VSWC_NEON_z := {
      
      x <- VSWC_NEON
      
      
      x_sd <- stats::sd(
        x
      )
      
      
      if (
        .N >= 3L &&
        is.finite(x_sd) &&
        x_sd > 0
      ) {
        
        (
          x -
            mean(
              x
            )
        ) /
          x_sd
        
      } else {
        
        rep(
          NA_real_,
          .N
        )
      }
    },
    by = .(
      index,
      verticalPosition
    )
  ]
  
  
  comparison[
    paired == TRUE,
    SoilWetnessFrac_SIPNET_z := {
      
      x <- SoilWetnessFrac_SIPNET
      
      
      x_sd <- stats::sd(
        x
      )
      
      
      if (
        .N >= 3L &&
        is.finite(x_sd) &&
        x_sd > 0
      ) {
        
        (
          x -
            mean(
              x
            )
        ) /
          x_sd
        
      } else {
        
        rep(
          NA_real_,
          .N
        )
      }
    },
    by = .(
      index,
      verticalPosition
    )
  ]
  
  
  # ==========================================================================
  # 15. Keep paired timestamps if requested
  # ==========================================================================
  
  if (
    isTRUE(
      keep_only_paired
    )
  ) {
    
    comparison <- comparison[
      paired == TRUE
    ]
  }
  
  
  # ==========================================================================
  # 16. Final column order
  # ==========================================================================
  
  preferred_columns <- c(
    "time",
    "index",
    "AmeriFlux_ID",
    "NEON_code",
    "NEON_name",
    "verticalPosition",
    "mean_zOffset",
    
    "VSWC_NEON",
    "SoilWetnessFrac_SIPNET",
    "SoilWater_SIPNET",
    
    "VSWC_NEON_z",
    "SoilWetnessFrac_SIPNET_z",
    
    "VSWC_NEON_temporal_sd",
    "VSWC_NEON_spatial_sd",
    
    "SoilWetnessFrac_SIPNET_sd",
    "SoilWater_SIPNET_sd",
    
    "n_30min",
    "mean_n_position",
    "n_SIPNET_ensemble",
    "n_paired_site_depth",
    
    "AirT_C",
    "Precip_mm_3h",
    "VPD_Pa",
    "RH",
    "WindSpeed",
    "PAR_umol_m2_s",
    
    "paired"
  )
  
  
  preferred_columns <- intersect(
    preferred_columns,
    names(
      comparison
    )
  )
  
  
  data.table::setcolorder(
    comparison,
    c(
      preferred_columns,
      setdiff(
        names(
          comparison
        ),
        preferred_columns
      )
    )
  )
  
  
  data.table::setorder(
    comparison,
    index,
    verticalPosition,
    time
  )
  
  
  # ==========================================================================
  # 17. Save
  # ==========================================================================
  
  if (!is.null(output_file)) {
    
    dir.create(
      dirname(
        output_file
      ),
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    data.table::fwrite(
      comparison,
      output_file
    )
    
    
    status_file <- sub(
      "\\.csv$",
      "_status.csv",
      output_file
    )
    
    
    data.table::fwrite(
      status,
      status_file
    )
  }
  
  
  # ==========================================================================
  # 18. Summary
  # ==========================================================================
  
  cat(
    "\n============================================\n",
    "SOIL MOISTURE VALIDATION TABLE\n",
    "============================================\n",
    "Requested indices:       ",
    length(
      target_indices
    ),
    "\n",
    
    "SIPNET indices found:   ",
    sum(
      status$sipnet_success
    ),
    "\n",
    
    "NEON site-depth pairs:  ",
    data.table::uniqueN(
      paste(
        comparison$index,
        comparison$verticalPosition
      )
    ),
    "\n",
    
    "Paired 3-hour rows:     ",
    sum(
      comparison$paired,
      na.rm = TRUE
    ),
    "\n",
    
    "Date range:             ",
    start_date,
    " to ",
    end_date,
    "\n",
    
    "============================================\n",
    sep = ""
  )
  
  
  return(
    comparison
  )
}
