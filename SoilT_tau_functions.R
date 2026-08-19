# =============================================================================
# Soil temperature tau calibration and validation functions
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})


# =============================================================================
# NEON soil-temperature depth metadata
# =============================================================================

NEON_NOMINAL_DEPTH_CM <- c(
  "501" = 2,
  "502" = 6,
  "503" = 16,
  "504" = 26
)


#' Parse timestamps in UTC
#'
#' Converts common NEON timestamp formats to `POSIXct` using UTC.
#'
#' @param x Character or character-coercible vector containing timestamps.
#'
#' @return A `POSIXct` vector in UTC.
#'
#' @md
#' @keywords internal
#' @author Yang Gu
parse_time_utc <- function(x) {
  
  x <- as.character(x)
  
  out <- suppressWarnings(
    as.POSIXct(
      x,
      tz = "UTC"
    )
  )
  
  if (all(is.na(out))) {
    
    out <- suppressWarnings(
      as.POSIXct(
        x,
        format = "%Y-%m-%dT%H:%M:%SZ",
        tz = "UTC"
      )
    )
  }
  
  if (all(is.na(out))) {
    
    out <- suppressWarnings(
      as.POSIXct(
        x,
        format = "%Y-%m-%d %H:%M:%S",
        tz = "UTC"
      )
    )
  }
  
  return(out)
}


#' Calculate a correlation safely
#'
#' Calculates Pearson correlation after removing non-finite pairs. Returns
#' `NA_real_` when fewer than three valid pairs are available or either input
#' has zero variance.
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#'
#' @return Numeric Pearson correlation or `NA_real_`.
#'
#' @md
#' @keywords internal
#' @author Yang Gu
safe_cor <- function(
    x,
    y
) {
  
  ok <- is.finite(x) &
    is.finite(y)
  
  if (
    sum(ok) < 3L ||
    stats::sd(x[ok]) == 0 ||
    stats::sd(y[ok]) == 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  return(
    stats::cor(
      x[ok],
      y[ok]
    )
  )
}


#' Calculate soil-temperature validation metrics
#'
#' Calculates predictive R2, bias, RMSE, MAE, and Pearson correlation for
#' paired observed and predicted values.
#'
#' @param observed Numeric vector of observations.
#' @param predicted Numeric vector of predictions.
#'
#' @return A one-row data.table containing `n`, `r2`, `bias`, `rmse`, `mae`,
#'   and `correlation`.
#'
#' @md
#' @keywords internal
#' @author Yang Gu
calc_metrics <- function(
    observed,
    predicted
) {
  
  ok <- is.finite(observed) &
    is.finite(predicted)
  
  observed <- observed[ok]
  
  predicted <- predicted[ok]
  
  n_pair <- length(
    observed
  )
  
  if (n_pair == 0L) {
    
    return(
      data.table::data.table(
        n = 0L,
        r2 = NA_real_,
        bias = NA_real_,
        rmse = NA_real_,
        mae = NA_real_,
        correlation = NA_real_
      )
    )
  }
  
  residual <- predicted -
    observed
  
  denominator <- sum(
    (
      observed -
        mean(observed)
    )^2
  )
  
  r2 <- if (
    n_pair >= 2L &&
    is.finite(denominator) &&
    denominator > 0
  ) {
    
    1 -
      sum(
        residual^2
      ) /
      denominator
    
  } else {
    
    NA_real_
  }
  
  return(
    data.table::data.table(
      n =
        n_pair,
      
      r2 =
        r2,
      
      bias =
        mean(
          residual
        ),
      
      rmse =
        sqrt(
          mean(
            residual^2
          )
        ),
      
      mae =
        mean(
          abs(
            residual
          )
        ),
      
      correlation =
        safe_cor(
          observed,
          predicted
        )
    )
  )
}


#' Assign meteorological season from date
#'
#' Converts calendar month to meteorological season using DJF, MAM, JJA,
#' and SON.
#'
#' @param date Date or date-coercible vector.
#'
#' @return Character vector containing `DJF`, `MAM`, `JJA`, or `SON`.
#'
#' @md
#' @keywords internal
#' @author Yang Gu
season_from_date <- function(
    date
) {
  
  month <- as.integer(
    format(
      date,
      "%m"
    )
  )
  
  season <- rep(
    NA_character_,
    length(month)
  )
  
  season[
    month %in% c(
      12L,
      1L,
      2L
    )
  ] <- "DJF"
  
  season[
    month %in%
      3L:5L
  ] <- "MAM"
  
  season[
    month %in%
      6L:8L
  ] <- "JJA"
  
  season[
    month %in%
      9L:11L
  ] <- "SON"
  
  return(
    season
  )
}


#' Read a SIPNET climate forcing file
#'
#' Reads a SIPNET `.clim` file and standardizes columns used for
#' soil-temperature calibration. Both 12-column and legacy 14-column climate
#' formats are supported.
#'
#' @param path Path to a SIPNET `.clim` file.
#'
#' @return A data.table containing standardized climate forcing and derived
#'   `timestamp` and `date` columns.
#'
#' @md
#' @export
#' @author Yang Gu
read_sipnet_clim <- function(
    path
) {
  
  if (!file.exists(path)) {
    
    stop(
      "SIPNET .clim file does not exist: ",
      path,
      call. = FALSE
    )
  }
  
  clim <- data.table::fread(
    path,
    header = FALSE,
    data.table = TRUE,
    showProgress = FALSE
  )
  
  names_v2 <- c(
    "year",
    "doy",
    "hour",
    "timestep_days",
    "air_temp_C",
    "soil_temp_C",
    "par",
    "precip",
    "vpd",
    "vpd_soil",
    "canopy_vp",
    "wind_speed"
  )
  
  names_v1 <- c(
    "grid_index",
    names_v2,
    "soil_wetness"
  )
  
  if (ncol(clim) == 12L) {
    
    data.table::setnames(
      clim,
      names_v2
    )
    
  } else if (ncol(clim) == 14L) {
    
    data.table::setnames(
      clim,
      names_v1
    )
    
  } else {
    
    stop(
      "Expected 12 or 14 columns in .clim; found ",
      ncol(clim),
      ".",
      call. = FALSE
    )
  }
  
  required_columns <- c(
    "year",
    "doy",
    "hour",
    "timestep_days",
    "air_temp_C",
    "soil_temp_C"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(clim)
  )
  
  if (length(missing_columns) > 0L) {
    
    stop(
      "Could not identify required .clim columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  clim[
    ,
    `:=`(
      year =
        as.integer(
          year
        ),
      
      doy =
        as.numeric(
          doy
        ),
      
      hour =
        as.numeric(
          hour
        ),
      
      timestep_days =
        as.numeric(
          timestep_days
        ),
      
      air_temp_C =
        as.numeric(
          air_temp_C
        ),
      
      soil_temp_C =
        as.numeric(
          soil_temp_C
        )
    )
  ]
  
  clim[
    ,
    timestamp :=
      as.POSIXct(
        sprintf(
          "%04d-01-01 00:00:00",
          year
        ),
        tz = "UTC"
      ) +
      (doy - 1) *
      86400 +
      hour *
      3600
  ]
  
  clim[
    ,
    date :=
      as.Date(
        timestamp
      )
  ]
  
  data.table::setorder(
    clim,
    timestamp
  )
  
  if (nrow(clim) == 0L) {
    
    stop(
      "No records were read from the .clim file.",
      call. = FALSE
    )
  }
  
  if (
    stats::median(
      clim$soil_temp_C,
      na.rm = TRUE
    ) < -100
  ) {
    
    stop(
      paste0(
        "The .clim SoilT median is below -100 C; ",
        "check temperature units."
      ),
      call. = FALSE
    )
  }
  
  return(
    clim
  )
}


#' Read NEON soil temperature for a specified vertical position
#'
#' Reads monthly NEON environmental measurement files for a single site,
#' selects one soil-temperature vertical position, applies the NEON final
#' quality flag when available, averages horizontal positions, and aggregates
#' observations to daily means.
#'
#' @param site Character NEON site code.
#' @param env_dir Directory containing monthly NEON environmental Rda files.
#' @param start_year First year to read.
#' @param end_year Last year to read.
#' @param vertical_position Character NEON soil-temperature vertical position.
#' @param nominal_depth_cm Nominal sensor depth in centimeters.
#' @param min_plots_per_timestamp Minimum number of horizontal positions
#'   required at a timestamp.
#' @param min_timestamps_per_day Minimum number of valid subdaily timestamps
#'   required to retain a daily mean.
#'
#' @return A daily data.table containing observed soil temperature and coverage
#'   diagnostics. Sensor metadata are attached as attributes.
#'
#' @md
#' @export
#' @author Yang Gu
read_neon_shallow_soil_temperature <- function(
    site,
    env_dir,
    start_year,
    end_year,
    vertical_position,
    nominal_depth_cm,
    min_plots_per_timestamp = 1L,
    min_timestamps_per_day = 12L
) {
  
  monthly <- list()
  
  counter <- 1L
  
  available_positions <- character()
  
  zoffset_values_m <- numeric()
  
  for (
    year_i in
    seq.int(
      as.integer(
        start_year
      ),
      as.integer(
        end_year
      )
    )
  ) {
    
    for (month_i in 1:12) {
      
      file_i <- file.path(
        env_dir,
        sprintf(
          "env-meas-%s-%04d-%02d.Rda",
          site,
          year_i,
          month_i
        )
      )
      
      if (!file.exists(file_i)) {
        
        next
      }
      
      env_i <- new.env(
        parent =
          emptyenv()
      )
      
      load(
        file_i,
        envir =
          env_i
      )
      
      if (!exists(
        "site_data",
        envir = env_i,
        inherits = FALSE
      )) {
        
        next
      }
      
      site_data <- get(
        "site_data",
        envir = env_i,
        inherits = FALSE
      )
      
      level_two <- tryCatch(
        getElement(
          site_data,
          2L
        ),
        error = function(e) {
          NULL
        }
      )
      
      if (is.null(level_two)) {
        
        next
      }
      
      soil_temperature_raw <- tryCatch(
        getElement(
          level_two,
          3L
        ),
        error = function(e) {
          NULL
        }
      )
      
      if (is.null(soil_temperature_raw)) {
        
        next
      }
      
      dt <- tryCatch(
        data.table::as.data.table(
          soil_temperature_raw
        ),
        error = function(e) {
          NULL
        }
      )
      
      if (
        is.null(dt) ||
        nrow(dt) == 0L
      ) {
        
        next
      }
      
      required_columns <- c(
        "startDateTime",
        "horizontalPosition",
        "verticalPosition",
        "soilTempMean"
      )
      
      if (
        !all(
          required_columns %in%
          names(dt)
        )
      ) {
        
        next
      }
      
      dt[
        ,
        verticalPosition :=
          as.character(
            verticalPosition
          )
      ]
      
      dt[
        ,
        horizontalPosition :=
          as.character(
            horizontalPosition
          )
      ]
      
      available_positions <- union(
        available_positions,
        unique(
          dt$verticalPosition
        )
      )
      
      dt <- dt[
        verticalPosition ==
          as.character(
            vertical_position
          )
      ]
      
      if (nrow(dt) == 0L) {
        
        next
      }
      
      dt[
        ,
        timestamp :=
          parse_time_utc(
            startDateTime
          )
      ]
      
      dt[
        ,
        soilTempMean :=
          suppressWarnings(
            as.numeric(
              soilTempMean
            )
          )
      ]
      
      if (
        "soilTempFinalQF" %in%
        names(dt)
      ) {
        
        dt[
          ,
          soilTempFinalQF_num :=
            suppressWarnings(
              as.numeric(
                soilTempFinalQF
              )
            )
        ]
        
        dt <- dt[
          is.na(
            soilTempFinalQF_num
          ) |
            soilTempFinalQF_num == 0
        ]
      }
      
      if (
        "zOffset" %in%
        names(dt)
      ) {
        
        dt[
          ,
          zOffset_m :=
            suppressWarnings(
              as.numeric(
                zOffset
              )
            )
        ]
        
        zoffset_values_m <- c(
          zoffset_values_m,
          
          dt[
            is.finite(
              zOffset_m
            ),
            zOffset_m
          ]
        )
        
      } else {
        
        dt[
          ,
          zOffset_m :=
            NA_real_
        ]
      }
      
      dt <- dt[
        !is.na(
          timestamp
        ) &
          is.finite(
            soilTempMean
          )
      ]
      
      if (nrow(dt) == 0L) {
        
        next
      }
      
      # -----------------------------------------------------------------------
      # Average replicate observations within timestamp x horizontal position
      # -----------------------------------------------------------------------
      
      by_plot <- dt[
        ,
        .(
          Tsoil_plot_C =
            mean(
              soilTempMean,
              na.rm = TRUE
            ),
          
          zOffset_plot_m = {
            
            z_use <- zOffset_m[
              is.finite(
                zOffset_m
              )
            ]
            
            if (length(z_use) > 0L) {
              
              stats::median(
                z_use
              )
              
            } else {
              
              NA_real_
            }
          }
        ),
        by = .(
          timestamp,
          horizontalPosition
        )
      ]
      
      # -----------------------------------------------------------------------
      # Average horizontal positions at each timestamp
      # -----------------------------------------------------------------------
      
      by_timestamp <- by_plot[
        ,
        .(
          Tsoil_shallow_C =
            mean(
              Tsoil_plot_C,
              na.rm = TRUE
            ),
          
          n_plots =
            data.table::uniqueN(
              horizontalPosition
            ),
          
          mean_zOffset_m = {
            
            z_use <- zOffset_plot_m[
              is.finite(
                zOffset_plot_m
              )
            ]
            
            if (length(z_use) > 0L) {
              
              mean(
                z_use
              )
              
            } else {
              
              NA_real_
            }
          }
        ),
        by =
          timestamp
      ]
      
      by_timestamp <- by_timestamp[
        n_plots >=
          as.integer(
            min_plots_per_timestamp
          )
      ]
      
      if (nrow(by_timestamp) == 0L) {
        
        next
      }
      
      # -----------------------------------------------------------------------
      # Aggregate to daily means
      # -----------------------------------------------------------------------
      
      daily <- by_timestamp[
        ,
        .(
          Tsoil_obs_C =
            mean(
              Tsoil_shallow_C,
              na.rm = TRUE
            ),
          
          n_tsoil_timestamps =
            data.table::uniqueN(
              timestamp
            ),
          
          mean_n_plots =
            mean(
              n_plots,
              na.rm = TRUE
            ),
          
          min_n_plots =
            min(
              n_plots,
              na.rm = TRUE
            ),
          
          max_n_plots =
            max(
              n_plots,
              na.rm = TRUE
            ),
          
          mean_zOffset_m = {
            
            z_use <- mean_zOffset_m[
              is.finite(
                mean_zOffset_m
              )
            ]
            
            if (length(z_use) > 0L) {
              
              mean(
                z_use
              )
              
            } else {
              
              NA_real_
            }
          }
        ),
        by = .(
          date =
            as.Date(
              timestamp
            )
        )
      ]
      
      daily <- daily[
        n_tsoil_timestamps >=
          as.integer(
            min_timestamps_per_day
          )
      ]
      
      if (nrow(daily) == 0L) {
        
        next
      }
      
      monthly[counter] <- list(
        daily
      )
      
      counter <- counter +
        1L
    }
  }
  
  if (length(monthly) == 0L) {
    
    stop(
      paste0(
        "No usable NEON SoilT found for ",
        site,
        " at VER=",
        vertical_position,
        "."
      ),
      call. = FALSE
    )
  }
  
  out <- data.table::rbindlist(
    monthly,
    use.names = TRUE,
    fill = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Protect against duplicated month-boundary dates
  # ---------------------------------------------------------------------------
  
  out <- out[
    ,
    .(
      Tsoil_obs_C =
        stats::weighted.mean(
          Tsoil_obs_C,
          w =
            n_tsoil_timestamps,
          na.rm = TRUE
        ),
      
      n_tsoil_timestamps =
        sum(
          n_tsoil_timestamps,
          na.rm = TRUE
        ),
      
      mean_n_plots =
        stats::weighted.mean(
          mean_n_plots,
          w =
            n_tsoil_timestamps,
          na.rm = TRUE
        ),
      
      min_n_plots =
        min(
          min_n_plots,
          na.rm = TRUE
        ),
      
      max_n_plots =
        max(
          max_n_plots,
          na.rm = TRUE
        ),
      
      mean_zOffset_m = {
        
        z_use <- mean_zOffset_m[
          is.finite(
            mean_zOffset_m
          )
        ]
        
        if (length(z_use) > 0L) {
          
          mean(
            z_use
          )
          
        } else {
          
          NA_real_
        }
      }
    ),
    by =
      date
  ]
  
  data.table::setorder(
    out,
    date
  )
  
  attr(
    out,
    "vertical_position"
  ) <- as.character(
    vertical_position
  )
  
  attr(
    out,
    "nominal_depth_cm"
  ) <- as.numeric(
    nominal_depth_cm
  )
  
  attr(
    out,
    "available_vertical_positions"
  ) <- sort(
    available_positions
  )
  
  attr(
    out,
    "zoffset_m_median"
  ) <- if (
    length(zoffset_values_m) > 0L
  ) {
    
    stats::median(
      zoffset_values_m,
      na.rm = TRUE
    )
    
  } else {
    
    NA_real_
  }
  
  message(
    "NEON SoilT days retained: ",
    nrow(out)
  )
  
  message(
    "Date range: ",
    min(out$date),
    " to ",
    max(out$date)
  )
  
  message(
    "Median timestamps/day: ",
    stats::median(
      out$n_tsoil_timestamps,
      na.rm = TRUE
    ),
    " / 48"
  )
  
  message(
    "Median plots/timestamp: ",
    round(
      stats::median(
        out$mean_n_plots,
        na.rm = TRUE
      ),
      2
    )
  )
  
  return(
    out
  )
}


#' Predict soil temperature using a causal exponential filter
#'
#' Applies a first-order thermal-memory model to air temperature to generate
#' a soil-temperature time series.
#'
#' \deqn{
#' T_{soil,t} =
#' T_{soil,t-1} +
#' \left(1 - \exp(-\Delta t / \tau)\right)
#' \left(T_{air,t} - T_{soil,t-1}\right)
#' }
#'
#' @param tair_C Numeric vector of air temperature in degrees Celsius.
#' @param timestep_days Numeric vector giving timestep duration in days.
#' @param tau_days Positive soil thermal time constant in days.
#' @param initial_soil_temp_C Optional initial soil temperature in degrees C.
#'
#' @return Numeric vector of predicted soil temperature in degrees Celsius.
#'
#' @md
#' @export
#' @author Yang Gu
causal_exponential_filter <- function(
    tair_C,
    timestep_days,
    tau_days,
    initial_soil_temp_C = NULL
) {
  
  if (
    length(tau_days) != 1L ||
    !is.finite(tau_days) ||
    tau_days <= 0
  ) {
    
    stop(
      "`tau_days` must be one positive finite value.",
      call. = FALSE
    )
  }
  
  n_timestep <- length(
    tair_C
  )
  
  if (n_timestep == 0L) {
    
    return(
      numeric()
    )
  }
  
  if (length(timestep_days) == 1L) {
    
    timestep_days <- rep(
      timestep_days,
      n_timestep
    )
  }
  
  if (
    length(timestep_days) !=
    n_timestep
  ) {
    
    stop(
      paste0(
        "`timestep_days` must have length 1 or the same length as ",
        "`tair_C`."
      ),
      call. = FALSE
    )
  }
  
  valid_timestep <- timestep_days[
    is.finite(
      timestep_days
    ) &
      timestep_days > 0
  ]
  
  if (length(valid_timestep) == 0L) {
    
    stop(
      "No positive finite timestep values were found.",
      call. = FALSE
    )
  }
  
  typical_timestep <- stats::median(
    valid_timestep
  )
  
  tsoil_C <- rep(
    NA_real_,
    n_timestep
  )
  
  first_valid_candidates <- which(
    is.finite(
      tair_C
    )
  )
  
  if (
    length(
      first_valid_candidates
    ) == 0L
  ) {
    
    return(
      tsoil_C
    )
  }
  
  first_valid <-
    first_valid_candidates[1L]
  
  tsoil_C[first_valid] <- if (
    is.null(
      initial_soil_temp_C
    )
  ) {
    
    tair_C[first_valid]
    
  } else {
    
    initial_value <- as.numeric(
      initial_soil_temp_C
    )[1L]
    
    if (
      is.finite(
        initial_value
      )
    ) {
      
      initial_value
      
    } else {
      
      tair_C[first_valid]
    }
  }
  
  if (
    first_valid <
    n_timestep
  ) {
    
    for (
      timestep_i in
      seq.int(
        first_valid + 1L,
        n_timestep
      )
    ) {
      
      previous_tsoil <- tsoil_C[
        timestep_i - 1L
      ]
      
      if (
        !is.finite(
          previous_tsoil
        )
      ) {
        
        previous_tsoil <- tair_C[
          timestep_i - 1L
        ]
      }
      
      if (
        !is.finite(
          tair_C[timestep_i]
        )
      ) {
        
        tsoil_C[timestep_i] <-
          previous_tsoil
        
        next
      }
      
      timestep_i_days <-
        timestep_days[
          timestep_i
        ]
      
      if (
        !is.finite(
          timestep_i_days
        ) ||
        timestep_i_days <= 0
      ) {
        
        timestep_i_days <-
          typical_timestep
      }
      
      alpha_i <-
        1 -
        exp(
          -timestep_i_days /
            tau_days
        )
      
      tsoil_C[timestep_i] <-
        previous_tsoil +
        alpha_i *
        (
          tair_C[timestep_i] -
            previous_tsoil
        )
    }
  }
  
  return(
    tsoil_C
  )
}


#' Predict daily soil temperature for a specified tau
#'
#' Runs the causal soil-temperature filter at the native SIPNET forcing
#' timestep and aggregates predicted soil temperature to daily means.
#'
#' @param clim A data.table containing `date`, `air_temp_C`, and
#'   `timestep_days`.
#' @param tau_days Positive soil thermal time constant in days.
#' @param offset_C Optional additive soil-temperature offset in degrees C.
#'
#' @return A data.table with `date` and `Tsoil_tau_C`.
#'
#' @md
#' @export
#' @author Yang Gu
predict_daily_tsoil <- function(
    clim,
    tau_days,
    offset_C = 0
) {
  
  required_columns <- c(
    "date",
    "air_temp_C",
    "timestep_days"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(clim)
  )
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "Missing required climate columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  predicted_subdaily <-
    causal_exponential_filter(
      tair_C =
        clim$air_temp_C,
      
      timestep_days =
        clim$timestep_days,
      
      tau_days =
        tau_days
    ) +
    offset_C
  
  prediction <- data.table::data.table(
    date =
      clim$date,
    
    Tsoil_tau_C =
      predicted_subdaily
  )[
    ,
    .(
      Tsoil_tau_C =
        mean(
          Tsoil_tau_C,
          na.rm = TRUE
        )
    ),
    by =
      date
  ]
  
  data.table::setorder(
    prediction,
    date
  )
  
  return(
    prediction
  )
}


#' Prepare daily soil-temperature data for tau calibration
#'
#' Aggregates SIPNET forcing to daily means, joins forcing with observed NEON
#' soil temperature, restricts data to requested observation years, removes
#' the initial forcing warm-up period, and assigns meteorological seasons.
#'
#' @param clim A data.table containing `date`, `air_temp_C`, and `soil_temp_C`.
#' @param obs_tsoil A data.table containing `date` and `Tsoil_obs_C`.
#' @param obs_start_year First observation year to retain.
#' @param obs_end_year Last observation year to retain.
#' @param warmup_days Number of forcing days used as initial thermal spin-up.
#'
#' @return A daily analysis data.table.
#'
#' @md
#' @export
#' @author Yang Gu
prepare_tau_analysis_data <- function(
    clim,
    obs_tsoil,
    obs_start_year,
    obs_end_year,
    warmup_days = 180L
) {
  
  required_clim_columns <- c(
    "date",
    "air_temp_C",
    "soil_temp_C"
  )
  
  required_obs_columns <- c(
    "date",
    "Tsoil_obs_C"
  )
  
  missing_clim <- setdiff(
    required_clim_columns,
    names(clim)
  )
  
  missing_obs <- setdiff(
    required_obs_columns,
    names(obs_tsoil)
  )
  
  if (length(missing_clim) > 0L) {
    
    stop(
      "Missing climate columns: ",
      paste(
        missing_clim,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (length(missing_obs) > 0L) {
    
    stop(
      "Missing observed SoilT columns: ",
      paste(
        missing_obs,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  model_daily <- clim[
    ,
    .(
      Tair_C =
        mean(
          air_temp_C,
          na.rm = TRUE
        ),
      
      Tsoil_current_C =
        mean(
          soil_temp_C,
          na.rm = TRUE
        )
    ),
    by =
      date
  ]
  
  analysis_data <- merge(
    model_daily,
    obs_tsoil,
    by = "date",
    all = FALSE
  )
  
  data.table::setorder(
    analysis_data,
    date
  )
  
  analysis_data[
    ,
    obs_year :=
      as.integer(
        format(
          date,
          "%Y"
        )
      )
  ]
  
  analysis_data <- analysis_data[
    obs_year >=
      as.integer(
        obs_start_year
      ) &
      obs_year <=
      as.integer(
        obs_end_year
      )
  ]
  
  forcing_warmup_end <-
    min(
      clim$date,
      na.rm = TRUE
    ) +
    as.integer(
      warmup_days
    )
  
  analysis_data <- analysis_data[
    date >=
      forcing_warmup_end
  ]
  
  if (
    nrow(
      analysis_data
    ) < 100L
  ) {
    
    stop(
      paste0(
        "Fewer than 100 overlapping daily SoilT observations remain ",
        "after observation-year and warm-up filtering."
      ),
      call. = FALSE
    )
  }
  
  analysis_data[
    ,
    season :=
      season_from_date(
        date
      )
  ]
  
  return(
    analysis_data
  )
}


#' Evaluate the soil-temperature likelihood for a candidate tau
#'
#' Calculates daily SoilT for a candidate thermal time constant and evaluates
#' a Gaussian likelihood against observed NEON soil temperature.
#'
#' @param clim Meteorological forcing used by `predict_daily_tsoil()`.
#' @param obs_daily Daily observed soil temperature.
#' @param train_dates Dates used to evaluate the likelihood.
#' @param tau_days Positive candidate tau in days.
#' @param fit_offset Logical indicating whether an additive offset is profiled.
#' @param min_fit_days Minimum number of paired daily observations required.
#'
#' @return A list containing `tau_days`, `offset_C`, `sigma_C`, `nll`, and `n`.
#'
#' @md
#' @export
#' @author Yang Gu
evaluate_tau_likelihood <- function(
    clim,
    obs_daily,
    train_dates,
    tau_days,
    fit_offset = FALSE,
    min_fit_days = 30L
) {
  
  train_dates <- unique(
    as.Date(
      train_dates
    )
  )
  
  obs_train <- obs_daily[
    date %in%
      train_dates
  ]
  
  prediction <- predict_daily_tsoil(
    clim =
      clim,
    
    tau_days =
      tau_days,
    
    offset_C =
      0
  )
  
  dat <- merge(
    obs_train,
    prediction,
    by = "date",
    all = FALSE
  )
  
  dat <- dat[
    is.finite(
      Tsoil_obs_C
    ) &
      is.finite(
        Tsoil_tau_C
      )
  ]
  
  if (
    nrow(dat) <
    as.integer(
      min_fit_days
    )
  ) {
    
    return(
      list(
        tau_days =
          tau_days,
        
        offset_C =
          NA_real_,
        
        sigma_C =
          NA_real_,
        
        nll =
          Inf,
        
        n =
          nrow(dat)
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Profile the additive temperature offset
  #
  # offset_MLE = mean(Tobs - Tpred)
  # ---------------------------------------------------------------------------
  
  offset_hat <- if (
    isTRUE(
      fit_offset
    )
  ) {
    
    mean(
      dat$Tsoil_obs_C -
        dat$Tsoil_tau_C,
      na.rm = TRUE
    )
    
  } else {
    
    0
  }
  
  residual <-
    dat$Tsoil_obs_C -
    (
      dat$Tsoil_tau_C +
        offset_hat
    )
  
  # ---------------------------------------------------------------------------
  # Profile Gaussian residual variance
  #
  # sigma^2_MLE = SSE / n
  # ---------------------------------------------------------------------------
  
  sigma_hat <- sqrt(
    mean(
      residual^2,
      na.rm = TRUE
    )
  )
  
  if (
    !is.finite(
      sigma_hat
    ) ||
    sigma_hat <= 0
  ) {
    
    return(
      list(
        tau_days =
          tau_days,
        
        offset_C =
          offset_hat,
        
        sigma_C =
          NA_real_,
        
        nll =
          Inf,
        
        n =
          nrow(dat)
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Gaussian negative log-likelihood
  # ---------------------------------------------------------------------------
  
  nll <- 0.5 *
    sum(
      log(
        2 *
          pi *
          sigma_hat^2
      ) +
        (
          residual /
            sigma_hat
        )^2
    )
  
  return(
    list(
      tau_days =
        tau_days,
      
      offset_C =
        offset_hat,
      
      sigma_C =
        sigma_hat,
      
      nll =
        nll,
      
      n =
        nrow(dat)
    )
  )
}


#' Estimate the soil thermal time constant by maximum likelihood
#'
#' Estimates tau by minimizing the profiled Gaussian negative log-likelihood
#' of observed NEON soil temperature. Optimization is performed on `log(tau)`
#' so the estimated thermal time constant remains positive.
#'
#' @param clim Meteorological forcing.
#' @param obs_daily Daily observed soil temperature.
#' @param train_dates Dates used for fitting.
#' @param tau_lower Positive lower tau bound in days.
#' @param tau_upper Positive upper tau bound in days.
#' @param fit_offset Logical indicating whether to profile an additive offset.
#' @param min_fit_days Minimum number of paired observations required.
#'
#' @return A list containing the MLE tau and likelihood diagnostics.
#'
#' @md
#' @export
#' @author Yang Gu
fit_tau_mle <- function(
    clim,
    obs_daily,
    train_dates,
    tau_lower = 0.5,
    tau_upper = 180,
    fit_offset = FALSE,
    min_fit_days = 30L
) {
  
  train_dates <- unique(
    as.Date(
      train_dates
    )
  )
  
  train_dates <- train_dates[
    !is.na(
      train_dates
    )
  ]
  
  if (
    length(train_dates) <
    as.integer(
      min_fit_days
    )
  ) {
    
    stop(
      "At least ",
      min_fit_days,
      " training dates are required for tau fitting.",
      call. = FALSE
    )
  }
  
  if (
    length(tau_lower) != 1L ||
    length(tau_upper) != 1L ||
    !is.finite(tau_lower) ||
    !is.finite(tau_upper) ||
    tau_lower <= 0 ||
    tau_upper <= tau_lower
  ) {
    
    stop(
      "`tau_lower` and `tau_upper` must define a positive increasing interval.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Optimize theta = log(tau)
  # ---------------------------------------------------------------------------
  
  objective_log_tau <- function(
    log_tau
  ) {
    
    tau_i <- exp(
      log_tau
    )
    
    likelihood_i <- evaluate_tau_likelihood(
      clim =
        clim,
      
      obs_daily =
        obs_daily,
      
      train_dates =
        train_dates,
      
      tau_days =
        tau_i,
      
      fit_offset =
        fit_offset,
      
      min_fit_days =
        min_fit_days
    )
    
    return(
      likelihood_i$nll
    )
  }
  
  fit <- stats::optimize(
    f =
      objective_log_tau,
    
    interval = c(
      log(
        tau_lower
      ),
      
      log(
        tau_upper
      )
    ),
    
    tol =
      1e-7
  )
  
  if (
    !is.finite(
      fit$objective
    )
  ) {
    
    stop(
      "Tau optimization did not identify a finite likelihood.",
      call. = FALSE
    )
  }
  
  tau_hat <- exp(
    fit$minimum
  )
  
  best <- evaluate_tau_likelihood(
    clim =
      clim,
    
    obs_daily =
      obs_daily,
    
    train_dates =
      train_dates,
    
    tau_days =
      tau_hat,
    
    fit_offset =
      fit_offset,
    
    min_fit_days =
      min_fit_days
  )
  
  return(
    list(
      tau_days =
        best$tau_days,
      
      offset_C =
        best$offset_C,
      
      sigma_C =
        best$sigma_C,
      
      nll =
        best$nll,
      
      n_train =
        best$n,
      
      convergence =
        0L,
      
      optimizer =
        "stats::optimize on log(tau) using profiled Gaussian NLL"
    )
  )
}


#' Evaluate the profile likelihood across a tau grid
#'
#' @param clim Meteorological forcing.
#' @param obs_daily Daily observed soil temperature.
#' @param train_dates Dates used in the likelihood calculation.
#' @param tau_values Positive vector of tau values in days.
#' @param fit_offset Logical indicating whether an additive offset is profiled.
#' @param min_fit_days Minimum number of paired observations required.
#'
#' @return A data.table containing likelihood values across tau.
#'
#' @md
#' @export
#' @author Yang Gu
profile_tau <- function(
    clim,
    obs_daily,
    train_dates,
    tau_values,
    fit_offset = FALSE,
    min_fit_days = 30L
) {
  
  tau_values <- as.numeric(
    tau_values
  )
  
  tau_values <- tau_values[
    is.finite(
      tau_values
    ) &
      tau_values > 0
  ]
  
  if (
    length(
      tau_values
    ) == 0L
  ) {
    
    stop(
      "`tau_values` must contain at least one positive finite value.",
      call. = FALSE
    )
  }
  
  profile_rows <- lapply(
    tau_values,
    
    function(
    tau_i
    ) {
      
      likelihood_i <- evaluate_tau_likelihood(
        clim =
          clim,
        
        obs_daily =
          obs_daily,
        
        train_dates =
          train_dates,
        
        tau_days =
          tau_i,
        
        fit_offset =
          fit_offset,
        
        min_fit_days =
          min_fit_days
      )
      
      data.table::data.table(
        tau_days =
          likelihood_i$tau_days,
        
        offset_C =
          likelihood_i$offset_C,
        
        sigma_C =
          likelihood_i$sigma_C,
        
        nll =
          likelihood_i$nll,
        
        n =
          likelihood_i$n
      )
    }
  )
  
  profile_table <- data.table::rbindlist(
    profile_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  data.table::setorder(
    profile_table,
    tau_days
  )
  
  return(
    profile_table
  )
}


#' Calculate an approximate profile-likelihood confidence interval for tau
#'
#' @param profile_table A table returned by `profile_tau()`.
#' @param confidence Confidence level.
#'
#' @return Numeric lower and upper confidence limits for tau.
#'
#' @md
#' @export
#' @author Yang Gu
profile_likelihood_ci <- function(
    profile_table,
    confidence = 0.95
) {
  
  required_columns <- c(
    "tau_days",
    "nll"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(
      profile_table
    )
  )
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "`profile_table` must contain `tau_days` and `nll`.",
      call. = FALSE
    )
  }
  
  if (
    length(confidence) != 1L ||
    !is.finite(confidence) ||
    confidence <= 0 ||
    confidence >= 1
  ) {
    
    stop(
      "`confidence` must be between 0 and 1.",
      call. = FALSE
    )
  }
  
  valid_profile <- profile_table[
    is.finite(
      tau_days
    ) &
      is.finite(
        nll
      )
  ]
  
  if (
    nrow(
      valid_profile
    ) == 0L
  ) {
    
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }
  
  nll_min <- min(
    valid_profile$nll
  )
  
  likelihood_threshold <-
    nll_min +
    0.5 *
    stats::qchisq(
      confidence,
      df = 1
    )
  
  inside_interval <- valid_profile[
    nll <=
      likelihood_threshold
  ]
  
  if (
    nrow(
      inside_interval
    ) == 0L
  ) {
    
    return(
      c(
        NA_real_,
        NA_real_
      )
    )
  }
  
  return(
    range(
      inside_interval$tau_days,
      na.rm = TRUE
    )
  )
}


#' Generate subdaily soil temperature using the final fitted tau
#'
#' @param clim SIPNET climate forcing.
#' @param final_fit Object returned by `fit_tau_mle()`.
#'
#' @return A data.table containing original and optimized subdaily SoilT.
#'
#' @md
#' @export
#' @author Yang Gu
make_corrected_subdaily <- function(
    clim,
    final_fit
) {
  
  required_columns <- c(
    "timestamp",
    "date",
    "air_temp_C",
    "soil_temp_C",
    "timestep_days"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(clim)
  )
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "Missing climate columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  tau_days <- getElement(
    final_fit,
    "tau_days"
  )
  
  offset_C <- getElement(
    final_fit,
    "offset_C"
  )
  
  if (
    is.null(tau_days) ||
    length(tau_days) != 1L ||
    !is.finite(tau_days) ||
    tau_days <= 0
  ) {
    
    stop(
      "`final_fit` does not contain a valid `tau_days` estimate.",
      call. = FALSE
    )
  }
  
  if (
    is.null(offset_C) ||
    length(offset_C) != 1L ||
    !is.finite(offset_C)
  ) {
    
    offset_C <- 0
  }
  
  optimized_tsoil <-
    causal_exponential_filter(
      tair_C =
        clim$air_temp_C,
      
      timestep_days =
        clim$timestep_days,
      
      tau_days =
        tau_days
    ) +
    offset_C
  
  return(
    data.table::data.table(
      timestamp =
        clim$timestamp,
      
      date =
        clim$date,
      
      air_temp_C =
        clim$air_temp_C,
      
      soil_temp_current_C =
        clim$soil_temp_C,
      
      soil_temp_tau_opt_C =
        optimized_tsoil
    )
  )
}


#' Summarize annual soil-temperature observation coverage
#'
#' @param analysis_data Data containing `date`, `obs_year`, and `season`.
#'
#' @return Annual and seasonal coverage table.
#'
#' @md
#' @export
#' @author Yang Gu
summarize_loyo_year_coverage <- function(
    analysis_data
) {
  
  required_columns <- c(
    "date",
    "obs_year",
    "season"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(
      analysis_data
    )
  )
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "Missing LOYO coverage columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  coverage <- analysis_data[
    ,
    .(
      n_days =
        data.table::uniqueN(
          date
        ),
      
      n_DJF =
        data.table::uniqueN(
          date[
            season == "DJF"
          ]
        ),
      
      n_MAM =
        data.table::uniqueN(
          date[
            season == "MAM"
          ]
        ),
      
      n_JJA =
        data.table::uniqueN(
          date[
            season == "JJA"
          ]
        ),
      
      n_SON =
        data.table::uniqueN(
          date[
            season == "SON"
          ]
        )
    ),
    by =
      obs_year
  ]
  
  data.table::setorder(
    coverage,
    obs_year
  )
  
  return(
    coverage
  )
}


#' Select years eligible for leave-one-year-out validation
#'
#' @param analysis_data Prepared daily analysis table.
#' @param min_days_per_year Minimum days required per year.
#' @param min_days_per_season Minimum days required per season.
#'
#' @return List containing coverage and eligible years.
#'
#' @md
#' @export
#' @author Yang Gu
select_loyo_years <- function(
    analysis_data,
    min_days_per_year = 240L,
    min_days_per_season = 30L
) {
  
  coverage <- summarize_loyo_year_coverage(
    analysis_data
  )
  
  coverage[
    ,
    eligible :=
      n_days >=
      as.integer(
        min_days_per_year
      ) &
      n_DJF >=
      as.integer(
        min_days_per_season
      ) &
      n_MAM >=
      as.integer(
        min_days_per_season
      ) &
      n_JJA >=
      as.integer(
        min_days_per_season
      ) &
      n_SON >=
      as.integer(
        min_days_per_season
      )
  ]
  
  return(
    list(
      coverage =
        coverage,
      
      eligible_years =
        coverage[
          eligible == TRUE,
          obs_year
        ]
    )
  )
}


#' Validate soil-temperature tau using leave-one-year-out cross-validation
#'
#' For each eligible year, tau is fitted without observations from that year
#' and then used to predict the held-out year.
#'
#' @param site Character NEON site identifier.
#' @param clim SIPNET meteorological forcing.
#' @param obs_tsoil Observed daily NEON SoilT.
#' @param analysis_data Prepared daily calibration data.
#' @param tau_lower Lower tau bound in days.
#' @param tau_upper Upper tau bound in days.
#' @param fit_offset Logical indicating whether an additive offset is profiled.
#' @param min_days_per_year Minimum annual coverage.
#' @param min_days_per_season Minimum seasonal coverage.
#' @param min_train_days Minimum training days per fold.
#'
#' @return A list containing LOYO predictions, tau estimates, and metrics.
#'
#' @md
#' @export
#' @author Yang Gu
run_loyo_tsoil_validation <- function(
    site,
    clim,
    obs_tsoil,
    analysis_data,
    tau_lower = 0.5,
    tau_upper = 180,
    fit_offset = FALSE,
    min_days_per_year = 240L,
    min_days_per_season = 30L,
    min_train_days = 180L
) {
  
  year_info <- select_loyo_years(
    analysis_data =
      analysis_data,
    
    min_days_per_year =
      min_days_per_year,
    
    min_days_per_season =
      min_days_per_season
  )
  
  eligible_years <- getElement(
    year_info,
    "eligible_years"
  )
  
  if (
    length(
      eligible_years
    ) < 2L
  ) {
    
    stop(
      "LOYO requires at least two eligible years. Eligible years: ",
      paste(
        eligible_years,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Fixed causal tau = 15 baseline
  # ---------------------------------------------------------------------------
  
  pred_tau15 <- predict_daily_tsoil(
    clim =
      clim,
    
    tau_days =
      15,
    
    offset_C =
      0
  )
  
  data.table::setnames(
    pred_tau15,
    "Tsoil_tau_C",
    "Tsoil_causal_tau15_C"
  )
  
  prediction_methods <- c(
    Current_SIPNET =
      "Tsoil_current_C",
    
    Causal_tau15 =
      "Tsoil_causal_tau15_C",
    
    Optimized_tau =
      "Tsoil_tau_opt_C"
  )
  
  tau_rows <- list()
  
  metric_rows <- list()
  
  validation_rows <- list()
  
  # ---------------------------------------------------------------------------
  # Leave one year out
  # ---------------------------------------------------------------------------
  
  for (
    held_out_year_i in
    eligible_years
  ) {
    
    train_dates <- analysis_data[
      obs_year !=
        held_out_year_i,
      unique(
        date
      )
    ]
    
    test_dates <- analysis_data[
      obs_year ==
        held_out_year_i,
      unique(
        date
      )
    ]
    
    if (
      length(train_dates) <
      as.integer(
        min_train_days
      )
    ) {
      
      warning(
        "Skipping held-out year ",
        held_out_year_i,
        ": only ",
        length(train_dates),
        " training days.",
        call. = FALSE
      )
      
      next
    }
    
    message(
      "LOYO ",
      site,
      " | held-out year = ",
      held_out_year_i,
      " | training days = ",
      length(train_dates),
      " | testing days = ",
      length(test_dates)
    )
    
    # -------------------------------------------------------------------------
    # Fit tau without the held-out year
    # -------------------------------------------------------------------------
    
    fold_fit <- fit_tau_mle(
      clim =
        clim,
      
      obs_daily =
        obs_tsoil,
      
      train_dates =
        train_dates,
      
      tau_lower =
        tau_lower,
      
      tau_upper =
        tau_upper,
      
      fit_offset =
        fit_offset,
      
      min_fit_days =
        min_train_days
    )
    
    # -------------------------------------------------------------------------
    # Predict with fold-specific optimized tau
    # -------------------------------------------------------------------------
    
    pred_optimized <- predict_daily_tsoil(
      clim =
        clim,
      
      tau_days =
        fold_fit$tau_days,
      
      offset_C =
        fold_fit$offset_C
    )
    
    data.table::setnames(
      pred_optimized,
      "Tsoil_tau_C",
      "Tsoil_tau_opt_C"
    )
    
    fold_data <- Reduce(
      function(
    x,
    y
      ) {
        
        merge(
          x,
          y,
          by = "date",
          all = FALSE
        )
      },
    
    list(
      analysis_data,
      pred_tau15,
      pred_optimized
    )
    )
    
    test_data <- fold_data[
      obs_year ==
        held_out_year_i
    ]
    
    if (
      nrow(
        test_data
      ) == 0L
    ) {
      
      next
    }
    
    test_data[
      ,
      held_out_year :=
        held_out_year_i
    ]
    
    # IMPORTANT:
    # Use single-bracket list assignment.
    # No nested double-bracket syntax is used anywhere.
    
    validation_rows[
      length(validation_rows) +
        1L
    ] <- list(
      test_data
    )
    
    tau_rows[
      length(tau_rows) +
        1L
    ] <- list(
      data.table::data.table(
        held_out_year =
          held_out_year_i,
        
        tau_days =
          fold_fit$tau_days,
        
        offset_C =
          fold_fit$offset_C,
        
        sigma_C =
          fold_fit$sigma_C,
        
        nll =
          fold_fit$nll,
        
        n_train =
          fold_fit$n_train
      )
    )
    
    # -------------------------------------------------------------------------
    # Held-out metrics
    # -------------------------------------------------------------------------
    
    for (
      method_name in
      names(
        prediction_methods
      )
    ) {
      
      prediction_column <- as.character(
        prediction_methods[
          method_name
        ]
      )
      
      predicted_values <- getElement(
        test_data,
        prediction_column
      )
      
      metrics <- calc_metrics(
        observed =
          test_data$Tsoil_obs_C,
        
        predicted =
          predicted_values
      )
      
      metrics[
        ,
        `:=`(
          held_out_year =
            held_out_year_i,
          
          method =
            method_name,
          
          fold_tau_days =
            fold_fit$tau_days
        )
      ]
      
      metric_rows[
        length(metric_rows) +
          1L
      ] <- list(
        metrics
      )
    }
  }
  
  if (
    length(
      validation_rows
    ) == 0L
  ) {
    
    stop(
      "No LOYO validation folds were completed.",
      call. = FALSE
    )
  }
  
  validation <- data.table::rbindlist(
    validation_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  tau_folds <- data.table::rbindlist(
    tau_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  fold_metrics <- data.table::rbindlist(
    metric_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Pooled out-of-sample metrics
  # ---------------------------------------------------------------------------
  
  pooled_rows <- list()
  
  for (
    method_name in
    names(
      prediction_methods
    )
  ) {
    
    prediction_column <- as.character(
      prediction_methods[
        method_name
      ]
    )
    
    predicted_values <- getElement(
      validation,
      prediction_column
    )
    
    metrics <- calc_metrics(
      observed =
        validation$Tsoil_obs_C,
      
      predicted =
        predicted_values
    )
    
    metrics[
      ,
      `:=`(
        method =
          method_name,
        
        validation_scheme =
          "leave_one_well_observed_year_out",
        
        n_folds =
          data.table::uniqueN(
            validation$held_out_year
          )
      )
    ]
    
    pooled_rows[
      length(pooled_rows) +
        1L
    ] <- list(
      metrics
    )
  }
  
  pooled_metrics <- data.table::rbindlist(
    pooled_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  return(
    list(
      year_coverage =
        getElement(
          year_info,
          "coverage"
        ),
      
      eligible_years =
        eligible_years,
      
      tau_folds =
        tau_folds,
      
      fold_metrics =
        fold_metrics,
      
      pooled_metrics =
        pooled_metrics,
      
      validation =
        validation,
      
      pred_tau15 =
        pred_tau15
    )
  )
}


#' Write tau calibration and validation outputs
#'
#' @param results Result list returned by `run_site_tau_workflow()`.
#' @param config Site-level configuration list.
#'
#' @return List containing summary table and output path.
#'
#' @md
#' @keywords internal
#' @author Yang Gu
write_tau_outputs <- function(
    results,
    config
) {
  
  output_dir <- getElement(
    config,
    "output_dir"
  )
  
  site <- getElement(
    config,
    "site"
  )
  
  obs_tsoil <- getElement(
    results,
    "obs_tsoil"
  )
  
  loyo <- getElement(
    results,
    "loyo"
  )
  
  final_fit <- getElement(
    results,
    "final_fit"
  )
  
  final_profile <- getElement(
    results,
    "final_profile"
  )
  
  final_profile_ci <- getElement(
    results,
    "final_profile_ci"
  )
  
  corrected_subdaily <- getElement(
    results,
    "corrected_subdaily"
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  data.table::fwrite(
    obs_tsoil,
    file.path(
      output_dir,
      paste0(
        site,
        "_NEON_shallow_tsoil_daily.csv"
      )
    )
  )
  
  data.table::fwrite(
    getElement(
      loyo,
      "year_coverage"
    ),
    file.path(
      output_dir,
      paste0(
        site,
        "_LOYO_year_coverage.csv"
      )
    )
  )
  
  data.table::fwrite(
    getElement(
      loyo,
      "validation"
    ),
    file.path(
      output_dir,
      paste0(
        site,
        "_soil_temperature_validation.csv"
      )
    )
  )
  
  data.table::fwrite(
    getElement(
      loyo,
      "pooled_metrics"
    ),
    file.path(
      output_dir,
      paste0(
        site,
        "_soil_temperature_metrics.csv"
      )
    )
  )
  
  data.table::fwrite(
    getElement(
      loyo,
      "fold_metrics"
    ),
    file.path(
      output_dir,
      paste0(
        site,
        "_soil_temperature_LOYO_fold_metrics.csv"
      )
    )
  )
  
  data.table::fwrite(
    getElement(
      loyo,
      "tau_folds"
    ),
    file.path(
      output_dir,
      paste0(
        site,
        "_tau_LOYO_folds.csv"
      )
    )
  )
  
  data.table::fwrite(
    final_profile,
    file.path(
      output_dir,
      paste0(
        site,
        "_tau_profile.csv"
      )
    )
  )
  
  data.table::fwrite(
    corrected_subdaily,
    file.path(
      output_dir,
      paste0(
        site,
        "_corrected_soil_temperature_subdaily.csv"
      )
    )
  )
  
  zoffset_m <- attr(
    obs_tsoil,
    "zoffset_m_median"
  )
  
  actual_depth_cm <- if (
    length(zoffset_m) == 1L &&
    is.finite(zoffset_m)
  ) {
    
    abs(zoffset_m) *
      100
    
  } else {
    
    NA_real_
  }
  
  pooled_metrics <- getElement(
    loyo,
    "pooled_metrics"
  )
  
  optimized_metrics <- pooled_metrics[
    method ==
      "Optimized_tau"
  ]
  
  tau_folds <- getElement(
    loyo,
    "tau_folds"
  )
  
  get_metric_value <- function(
    column_name
  ) {
    
    if (
      nrow(optimized_metrics) == 0L ||
      !column_name %in%
      names(
        optimized_metrics
      )
    ) {
      
      return(
        NA_real_
      )
    }
    
    value <- getElement(
      optimized_metrics,
      column_name
    )
    
    return(
      as.numeric(
        value[1L]
      )
    )
  }
  
  summary_table <- data.table::data.table(
    index =
      getElement(
        config,
        "index"
      ),
    
    AmeriFlux_ID =
      getElement(
        config,
        "AmeriFlux_ID"
      ),
    
    NEON_code =
      site,
    
    vertical_position =
      getElement(
        config,
        "vertical_position"
      ),
    
    nominal_depth_cm =
      getElement(
        config,
        "nominal_depth_cm"
      ),
    
    observed_zoffset_m_median =
      if (
        length(zoffset_m) ==
        1L
      ) {
        
        zoffset_m
        
      } else {
        
        NA_real_
      },
    
    observed_actual_depth_cm_median =
      actual_depth_cm,
    
    validation_scheme =
      "leave_one_well_observed_year_out",
    
    n_cv_years =
      nrow(
        tau_folds
      ),
    
    cv_tau_mean_days =
      mean(
        tau_folds$tau_days,
        na.rm = TRUE
      ),
    
    cv_tau_median_days =
      stats::median(
        tau_folds$tau_days,
        na.rm = TRUE
      ),
    
    cv_tau_sd_days =
      if (
        nrow(
          tau_folds
        ) > 1L
      ) {
        
        stats::sd(
          tau_folds$tau_days,
          na.rm = TRUE
        )
        
      } else {
        
        NA_real_
      },
    
    tau_days =
      getElement(
        final_fit,
        "tau_days"
      ),
    
    tau_ci_low_profile_approx =
      final_profile_ci[1L],
    
    tau_ci_high_profile_approx =
      final_profile_ci[2L],
    
    offset_C =
      getElement(
        final_fit,
        "offset_C"
      ),
    
    residual_sigma_C =
      getElement(
        final_fit,
        "sigma_C"
      ),
    
    final_nll =
      getElement(
        final_fit,
        "nll"
      ),
    
    LOYO_R2 =
      get_metric_value(
        "r2"
      ),
    
    LOYO_RMSE_C =
      get_metric_value(
        "rmse"
      ),
    
    LOYO_MAE_C =
      get_metric_value(
        "mae"
      ),
    
    LOYO_bias_C =
      get_metric_value(
        "bias"
      ),
    
    LOYO_correlation =
      get_metric_value(
        "correlation"
      )
  )
  
  summary_file <- file.path(
    output_dir,
    sprintf(
      "%s_index%d_VER%s_summary.csv",
      site,
      as.integer(
        getElement(
          config,
          "index"
        )
      ),
      as.character(
        getElement(
          config,
          "vertical_position"
        )
      )
    )
  )
  
  data.table::fwrite(
    summary_table,
    summary_file
  )
  
  return(
    list(
      summary =
        summary_table,
      
      summary_file =
        summary_file
    )
  )
}


#' Run the complete tau workflow for one NEON site and depth
#'
#' Reads forcing and observations, performs LOYO validation, fits final MLE
#' tau, calculates profile likelihood, generates optimized SoilT, and writes
#' outputs.
#'
#' @param config Site-depth configuration list.
#'
#' @return Invisibly returns complete result list.
#'
#' @md
#' @keywords internal
#' @author Yang Gu
run_site_tau_workflow <- function(
    config
) {
  
  output_dir <- getElement(
    config,
    "output_dir"
  )
  
  site <- getElement(
    config,
    "site"
  )
  
  vertical_position <- getElement(
    config,
    "vertical_position"
  )
  
  nominal_depth_cm <- getElement(
    config,
    "nominal_depth_cm"
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  message(
    "Site: ",
    site
  )
  
  message(
    "NEON target: VER=",
    vertical_position,
    " (nominal ~",
    nominal_depth_cm,
    " cm)"
  )
  
  message(
    "Reading SIPNET forcing: ",
    getElement(
      config,
      "clim_file"
    )
  )
  
  clim <- read_sipnet_clim(
    getElement(
      config,
      "clim_file"
    )
  )
  
  message(
    "Reading NEON SoilT from: ",
    getElement(
      config,
      "neon_env_dir"
    )
  )
  
  obs_tsoil <- read_neon_shallow_soil_temperature(
    site =
      site,
    
    env_dir =
      getElement(
        config,
        "neon_env_dir"
      ),
    
    start_year =
      getElement(
        config,
        "obs_start_year"
      ),
    
    end_year =
      getElement(
        config,
        "obs_end_year"
      ),
    
    vertical_position =
      vertical_position,
    
    nominal_depth_cm =
      nominal_depth_cm,
    
    min_plots_per_timestamp =
      getElement(
        config,
        "min_plots_per_timestamp"
      ),
    
    min_timestamps_per_day =
      getElement(
        config,
        "min_tsoil_timestamps_per_day"
      )
  )
  
  message(
    "Available NEON verticalPosition values: ",
    paste(
      attr(
        obs_tsoil,
        "available_vertical_positions"
      ),
      collapse = ", "
    )
  )
  
  zoffset_m <- attr(
    obs_tsoil,
    "zoffset_m_median"
  )
  
  if (
    length(zoffset_m) == 1L &&
    is.finite(zoffset_m)
  ) {
    
    message(
      sprintf(
        "Median zOffset: %.3f m (%.1f cm below surface).",
        zoffset_m,
        abs(zoffset_m) * 100
      )
    )
  }
  
  analysis_data <- prepare_tau_analysis_data(
    clim =
      clim,
    
    obs_tsoil =
      obs_tsoil,
    
    obs_start_year =
      getElement(
        config,
        "obs_start_year"
      ),
    
    obs_end_year =
      getElement(
        config,
        "obs_end_year"
      ),
    
    warmup_days =
      getElement(
        config,
        "warmup_days"
      )
  )
  
  # ---------------------------------------------------------------------------
  # LOYO validation
  # ---------------------------------------------------------------------------
  
  loyo <- run_loyo_tsoil_validation(
    site =
      site,
    
    clim =
      clim,
    
    obs_tsoil =
      obs_tsoil,
    
    analysis_data =
      analysis_data,
    
    tau_lower =
      getElement(
        config,
        "tau_lower_days"
      ),
    
    tau_upper =
      getElement(
        config,
        "tau_upper_days"
      ),
    
    fit_offset =
      getElement(
        config,
        "fit_temperature_offset"
      ),
    
    min_days_per_year =
      getElement(
        config,
        "min_loyo_days_per_year"
      ),
    
    min_days_per_season =
      getElement(
        config,
        "min_loyo_days_per_season"
      ),
    
    min_train_days =
      getElement(
        config,
        "min_train_tsoil_days"
      )
  )
  
  # ---------------------------------------------------------------------------
  # Final MLE using all available observation dates
  # ---------------------------------------------------------------------------
  
  all_dates <- analysis_data[
    ,
    unique(
      date
    )
  ]
  
  final_fit <- fit_tau_mle(
    clim =
      clim,
    
    obs_daily =
      obs_tsoil,
    
    train_dates =
      all_dates,
    
    tau_lower =
      getElement(
        config,
        "tau_lower_days"
      ),
    
    tau_upper =
      getElement(
        config,
        "tau_upper_days"
      ),
    
    fit_offset =
      getElement(
        config,
        "fit_temperature_offset"
      ),
    
    min_fit_days =
      getElement(
        config,
        "min_train_tsoil_days"
      )
  )
  
  # ---------------------------------------------------------------------------
  # Profile likelihood
  # ---------------------------------------------------------------------------
  
  tau_grid <- exp(
    seq(
      log(
        getElement(
          config,
          "tau_lower_days"
        )
      ),
      
      log(
        getElement(
          config,
          "tau_upper_days"
        )
      ),
      
      length.out =
        160L
    )
  )
  
  final_profile <- profile_tau(
    clim =
      clim,
    
    obs_daily =
      obs_tsoil,
    
    train_dates =
      all_dates,
    
    tau_values =
      tau_grid,
    
    fit_offset =
      getElement(
        config,
        "fit_temperature_offset"
      ),
    
    min_fit_days =
      getElement(
        config,
        "min_train_tsoil_days"
      )
  )
  
  final_profile_ci <- profile_likelihood_ci(
    final_profile,
    confidence = 0.95
  )
  
  corrected_subdaily <- make_corrected_subdaily(
    clim =
      clim,
    
    final_fit =
      final_fit
  )
  
  results <- list(
    config =
      config,
    
    clim =
      clim,
    
    obs_tsoil =
      obs_tsoil,
    
    analysis_data =
      analysis_data,
    
    loyo =
      loyo,
    
    final_fit =
      final_fit,
    
    final_profile =
      final_profile,
    
    final_profile_ci =
      final_profile_ci,
    
    corrected_subdaily =
      corrected_subdaily
  )
  
  output <- write_tau_outputs(
    results =
      results,
    
    config =
      config
  )
  
  results$summary <- getElement(
    output,
    "summary"
  )
  
  results$summary_file <- getElement(
    output,
    "summary_file"
  )
  
  tau_folds <- getElement(
    loyo,
    "tau_folds"
  )
  
  message(
    "Completed site-level tau workflow."
  )
  
  message(
    "LOYO held-out years: ",
    paste(
      tau_folds$held_out_year,
      collapse = ", "
    )
  )
  
  message(
    sprintf(
      "LOYO tau median: %.2f days",
      stats::median(
        tau_folds$tau_days,
        na.rm = TRUE
      )
    )
  )
  
  message(
    sprintf(
      "Final all-years tau: %.2f days",
      getElement(
        final_fit,
        "tau_days"
      )
    )
  )
  
  message(
    "Outputs: ",
    normalizePath(
      output_dir,
      mustWork = FALSE
    )
  )
  
  return(
    invisible(
      results
    )
  )
}


#' Estimate soil-temperature tau for one lookup index and depth
#'
#' Resolves a model index to a NEON site, identifies the SIPNET forcing file,
#' and runs MLE fitting and LOYO validation for one NEON soil depth.
#'
#' @param index Integer model-site index.
#' @param lookup Site lookup table containing `index` and `NEON_code`.
#' @param vertical_position Character NEON vertical position.
#' @param nominal_depth_cm Optional nominal depth in centimeters.
#' @param obs_start_year First observation year.
#' @param obs_end_year Last observation year.
#' @param clim_root Root directory containing SIPNET forcing.
#' @param clim_basename SIPNET climate file name.
#' @param neon_env_dir Directory containing NEON environmental Rda files.
#' @param output_root Output directory.
#' @param min_plots_per_timestamp Minimum horizontal plots per timestamp.
#' @param min_tsoil_timestamps_per_day Minimum observations per day.
#' @param tau_lower_days Lower MLE tau bound.
#' @param tau_upper_days Upper MLE tau bound.
#' @param fit_temperature_offset Logical indicating whether to fit offset.
#' @param warmup_days Initial model warm-up duration.
#' @param min_loyo_days_per_year Minimum annual LOYO coverage.
#' @param min_loyo_days_per_season Minimum seasonal LOYO coverage.
#' @param min_train_tsoil_days Minimum training observations.
#'
#' @return Invisibly returns summary, best tau, and complete result object.
#'
#' @md
#' @export
#' @author Yang Gu
estimate_tau_by_index <- function(
    index,
    lookup,
    vertical_position = "501",
    nominal_depth_cm = NULL,
    obs_start_year = 2017L,
    obs_end_year = 2021L,
    clim_root =
      "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024",
    clim_basename =
      "ERA5.1.2012-01-01.2024-12-31.clim",
    neon_env_dir =
      "/projectnb/dietzelab/jzobitz/02-NEON-sites/env-data",
    output_root =
      "/projectnb/dietzelab/guYANG/soilparam/tau_by_index",
    min_plots_per_timestamp = 1L,
    min_tsoil_timestamps_per_day = 12L,
    tau_lower_days = 0.5,
    tau_upper_days = 180,
    fit_temperature_offset = FALSE,
    warmup_days = 180L,
    min_loyo_days_per_year = 120L,
    min_loyo_days_per_season = 10L,
    min_train_tsoil_days = 60L
) {
  
  target_index <- suppressWarnings(
    as.integer(
      index
    )
  )
  
  if (
    length(target_index) != 1L ||
    is.na(target_index)
  ) {
    
    stop(
      "`index` must be one valid integer.",
      call. = FALSE
    )
  }
  
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
  
  required_lookup_columns <- c(
    "index",
    "NEON_code"
  )
  
  missing_lookup_columns <- setdiff(
    required_lookup_columns,
    names(
      lookup_dt
    )
  )
  
  if (
    length(
      missing_lookup_columns
    ) > 0L
  ) {
    
    stop(
      "lookup is missing: ",
      paste(
        missing_lookup_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  site_info <- lookup_dt[
    index ==
      target_index
  ]
  
  if (
    nrow(
      site_info
    ) == 0L
  ) {
    
    stop(
      "Index ",
      target_index,
      " was not found in lookup.",
      call. = FALSE
    )
  }
  
  if (
    nrow(
      site_info
    ) > 1L
  ) {
    
    stop(
      "Index ",
      target_index,
      " has multiple lookup rows.",
      call. = FALSE
    )
  }
  
  neon_code <- as.character(
    getElement(
      site_info,
      "NEON_code"
    )[1L]
  )
  
  ameriflux_id <- if (
    "AmeriFlux_ID" %in%
    names(
      site_info
    )
  ) {
    
    as.character(
      getElement(
        site_info,
        "AmeriFlux_ID"
      )[1L]
    )
    
  } else {
    
    NA_character_
  }
  
  if (
    is.na(neon_code) ||
    neon_code == ""
  ) {
    
    stop(
      "Index ",
      target_index,
      " has no valid NEON_code.",
      call. = FALSE
    )
  }
  
  vertical_position <- as.character(
    vertical_position
  )[1L]
  
  if (
    is.null(
      nominal_depth_cm
    )
  ) {
    
    if (
      !vertical_position %in%
      names(
        NEON_NOMINAL_DEPTH_CM
      )
    ) {
      
      stop(
        "No nominal depth is defined for vertical_position = ",
        vertical_position,
        ".",
        call. = FALSE
      )
    }
    
    nominal_depth_cm <- unname(
      NEON_NOMINAL_DEPTH_CM[
        vertical_position
      ]
    )
  }
  
  nominal_depth_cm <- as.numeric(
    nominal_depth_cm
  )[1L]
  
  if (
    !is.finite(
      nominal_depth_cm
    )
  ) {
    
    stop(
      "`nominal_depth_cm` must be finite.",
      call. = FALSE
    )
  }
  
  if (
    as.integer(
      obs_end_year
    ) <
    as.integer(
      obs_start_year
    )
  ) {
    
    stop(
      "`obs_end_year` must be >= `obs_start_year`.",
      call. = FALSE
    )
  }
  
  if (
    !is.finite(
      tau_lower_days
    ) ||
    !is.finite(
      tau_upper_days
    ) ||
    tau_lower_days <= 0 ||
    tau_upper_days <=
    tau_lower_days
  ) {
    
    stop(
      "Invalid tau bounds.",
      call. = FALSE
    )
  }
  
  clim_file <- file.path(
    clim_root,
    sprintf(
      "ERA5_%d_1",
      target_index
    ),
    clim_basename
  )
  
  if (
    !file.exists(
      clim_file
    )
  ) {
    
    stop(
      "CLIM file not found:\n",
      clim_file,
      call. = FALSE
    )
  }
  
  output_dir <- file.path(
    output_root,
    sprintf(
      "%s_index%d_VER%s",
      neon_code,
      target_index,
      vertical_position
    )
  )
  
  config <- list(
    index =
      target_index,
    
    AmeriFlux_ID =
      ameriflux_id,
    
    site =
      neon_code,
    
    obs_start_year =
      as.integer(
        obs_start_year
      ),
    
    obs_end_year =
      as.integer(
        obs_end_year
      ),
    
    clim_file =
      clim_file,
    
    neon_env_dir =
      neon_env_dir,
    
    vertical_position =
      vertical_position,
    
    nominal_depth_cm =
      nominal_depth_cm,
    
    min_plots_per_timestamp =
      as.integer(
        min_plots_per_timestamp
      ),
    
    min_tsoil_timestamps_per_day =
      as.integer(
        min_tsoil_timestamps_per_day
      ),
    
    tau_lower_days =
      as.numeric(
        tau_lower_days
      ),
    
    tau_upper_days =
      as.numeric(
        tau_upper_days
      ),
    
    fit_temperature_offset =
      isTRUE(
        fit_temperature_offset
      ),
    
    warmup_days =
      as.integer(
        warmup_days
      ),
    
    min_loyo_days_per_year =
      as.integer(
        min_loyo_days_per_year
      ),
    
    min_loyo_days_per_season =
      as.integer(
        min_loyo_days_per_season
      ),
    
    min_train_tsoil_days =
      as.integer(
        min_train_tsoil_days
      ),
    
    output_dir =
      output_dir
  )
  
  cat(
    "\n",
    "========================================\n",
    "SIPNET SoilT tau calibration\n",
    "========================================\n",
    "Index:       ", target_index, "\n",
    "NEON site:   ", neon_code, "\n",
    "AmeriFlux:   ", ameriflux_id, "\n",
    "Target VER:  ", vertical_position, "\n",
    "Depth:       ~", nominal_depth_cm, " cm\n",
    "Years:       ", obs_start_year, "-", obs_end_year, "\n",
    "========================================\n\n",
    sep = ""
  )
  
  results <- run_site_tau_workflow(
    config
  )
  
  summary_table <- getElement(
    results,
    "summary"
  )
  
  final_fit <- getElement(
    results,
    "final_fit"
  )
  
  print(
    summary_table
  )
  
  return(
    invisible(
      list(
        summary =
          summary_table,
        
        best_tau_days =
          getElement(
            final_fit,
            "tau_days"
          ),
        
        results =
          results
      )
    )
  )
}


#' Estimate tau for all lookup indices and NEON soil depths in parallel
#'
#' Runs `estimate_tau_by_index()` for every valid index in `lookup` and every
#' requested NEON soil-temperature depth.
#'
#' Parallelization occurs across model indices. Depths are processed
#' sequentially within each site.
#'
#' @param lookup Site lookup table.
#' @param workers Number of parallel workers.
#' @param vertical_positions NEON vertical positions to process.
#' @param depth_map Named vector mapping vertical position to nominal depth.
#' @param combined_output_file Output path for summary table.
#' @param status_output_file Output path for processing status.
#' @param ... Additional arguments passed to `estimate_tau_by_index()`.
#'
#' @return Invisibly returns combined summary and status tables.
#'
#' @md
#' @export
#' @author Yang Gu
run_all_indices_depths_parallel <- function(
    lookup,
    workers = 16L,
    vertical_positions = c(
      "501",
      "502",
      "503",
      "504"
    ),
    depth_map =
      NEON_NOMINAL_DEPTH_CM,
    combined_output_file =
      "/projectnb/dietzelab/guYANG/soilparam/tau_by_index/all_indices_all_depths_summary.csv",
    status_output_file =
      "/projectnb/dietzelab/guYANG/soilparam/tau_by_index/all_indices_all_depths_status.csv",
    ...
) {
  
  if (
    !requireNamespace(
      "future",
      quietly = TRUE
    )
  ) {
    
    stop(
      "Package 'future' is required.",
      call. = FALSE
    )
  }
  
  if (
    !requireNamespace(
      "future.apply",
      quietly = TRUE
    )
  ) {
    
    stop(
      "Package 'future.apply' is required.",
      call. = FALSE
    )
  }
  
  if (
    !requireNamespace(
      "data.table",
      quietly = TRUE
    )
  ) {
    
    stop(
      "Package 'data.table' is required.",
      call. = FALSE
    )
  }
  
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
  
  required_columns <- c(
    "index",
    "NEON_code"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(
      lookup_dt
    )
  )
  
  if (
    length(
      missing_columns
    ) > 0L
  ) {
    
    stop(
      "lookup is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  vertical_positions <- as.character(
    vertical_positions
  )
  
  if (
    !all(
      vertical_positions %in%
      names(
        depth_map
      )
    )
  ) {
    
    unknown_positions <- setdiff(
      vertical_positions,
      names(
        depth_map
      )
    )
    
    stop(
      "No nominal depth is defined for: ",
      paste(
        unknown_positions,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  index_values <- getElement(
    lookup_dt,
    "index"
  )
  
  indices <- sort(
    unique(
      as.integer(
        index_values
      )
    )
  )
  
  indices <- indices[
    is.finite(
      indices
    )
  ]
  
  if (
    length(
      indices
    ) == 0L
  ) {
    
    stop(
      "No valid lookup indices were found.",
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Capture additional arguments once.
  #
  # This avoids using `...` directly inside the parallel worker.
  # ---------------------------------------------------------------------------
  
  dots <- list(...)
  
  reserved_dot_names <- intersect(
    names(
      dots
    ),
    c(
      "index",
      "lookup",
      "vertical_position",
      "nominal_depth_cm"
    )
  )
  
  if (
    length(
      reserved_dot_names
    ) > 0L
  ) {
    
    stop(
      "Do not pass these arguments through `...`: ",
      paste(
        reserved_dot_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # One index per parallel job
  # ---------------------------------------------------------------------------
  
  run_one_index <- function(
    index_i
  ) {
    
    result_rows <- list()
    
    status_rows <- list()
    
    site_rows <- lookup_dt[
      index ==
        index_i
    ]
    
    if (
      nrow(
        site_rows
      ) != 1L
    ) {
      
      return(
        list(
          summary =
            NULL,
          
          status =
            data.table::data.table(
              index =
                index_i,
              
              NEON_code =
                NA_character_,
              
              vertical_position =
                NA_character_,
              
              nominal_depth_cm =
                NA_real_,
              
              status =
                "lookup_error",
              
              message =
                paste0(
                  "Expected one lookup row; found ",
                  nrow(
                    site_rows
                  )
                )
            )
        )
      )
    }
    
    neon_code <- as.character(
      getElement(
        site_rows,
        "NEON_code"
      )[1L]
    )
    
    for (
      vertical_position_i in
      vertical_positions
    ) {
      
      depth_i_cm <- as.numeric(
        depth_map[
          vertical_position_i
        ]
      )[1L]
      
      message(
        "\n----------------------------------------",
        "\nIndex: ",
        index_i,
        " | Site: ",
        neon_code,
        " | VER: ",
        vertical_position_i,
        " | Depth: ~",
        depth_i_cm,
        " cm",
        "\n----------------------------------------"
      )
      
      fit_args <- c(
        list(
          index =
            index_i,
          
          lookup =
            lookup_dt,
          
          vertical_position =
            vertical_position_i,
          
          nominal_depth_cm =
            depth_i_cm
        ),
        dots
      )
      
      fit_i <- tryCatch(
        do.call(
          estimate_tau_by_index,
          fit_args
        ),
        error = function(e) {
          e
        }
      )
      
      if (
        inherits(
          fit_i,
          "error"
        )
      ) {
        
        error_message <- conditionMessage(
          fit_i
        )
        
        no_depth_data <- grepl(
          paste0(
            "No usable NEON SoilT|",
            "No daily NEON SoilT|",
            "No usable SoilT|",
            "No NEON SoilT"
          ),
          error_message,
          ignore.case = TRUE
        )
        
        status_rows[
          length(status_rows) +
            1L
        ] <- list(
          data.table::data.table(
            index =
              index_i,
            
            NEON_code =
              neon_code,
            
            vertical_position =
              vertical_position_i,
            
            nominal_depth_cm =
              depth_i_cm,
            
            status =
              if (
                no_depth_data
              ) {
                
                "skipped_no_depth_data"
                
              } else {
                
                "failed"
              },
            
            message =
              error_message
          )
        )
        
        message(
          if (
            no_depth_data
          ) {
            
            "SKIP"
            
          } else {
            
            "FAILED"
          },
          " | index=",
          index_i,
          " | VER=",
          vertical_position_i,
          " | ",
          error_message
        )
        
        next
      }
      
      fit_summary <- getElement(
        fit_i,
        "summary"
      )
      
      best_tau_days <- getElement(
        fit_i,
        "best_tau_days"
      )
      
      result_i <- data.table::as.data.table(
        data.table::copy(
          fit_summary
        )
      )
      
      result_i[
        ,
        `:=`(
          vertical_position =
            vertical_position_i,
          
          nominal_depth_cm =
            depth_i_cm
        )
      ]
      
      result_rows[
        length(result_rows) +
          1L
      ] <- list(
        result_i
      )
      
      status_rows[
        length(status_rows) +
          1L
      ] <- list(
        data.table::data.table(
          index =
            index_i,
          
          NEON_code =
            neon_code,
          
          vertical_position =
            vertical_position_i,
          
          nominal_depth_cm =
            depth_i_cm,
          
          status =
            "success",
          
          message =
            NA_character_
        )
      )
      
      message(
        "SUCCESS | index=",
        index_i,
        " | VER=",
        vertical_position_i,
        " | tau=",
        round(
          best_tau_days,
          3
        ),
        " days"
      )
    }
    
    index_summary <- if (
      length(
        result_rows
      ) > 0L
    ) {
      
      data.table::rbindlist(
        result_rows,
        use.names = TRUE,
        fill = TRUE
      )
      
    } else {
      
      NULL
    }
    
    index_status <- if (
      length(
        status_rows
      ) > 0L
    ) {
      
      data.table::rbindlist(
        status_rows,
        use.names = TRUE,
        fill = TRUE
      )
      
    } else {
      
      data.table::data.table(
        index =
          index_i,
        
        NEON_code =
          neon_code,
        
        vertical_position =
          NA_character_,
        
        nominal_depth_cm =
          NA_real_,
        
        status =
          "no_result",
        
        message =
          NA_character_
      )
    }
    
    return(
      list(
        summary =
          index_summary,
        
        status =
          index_status
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Parallel execution across indices
  # ---------------------------------------------------------------------------
  
  future::plan(
    future::multisession,
    workers =
      as.integer(
        workers
      )
  )
  
  message(
    "Starting parallel tau calibration..."
  )
  
  message(
    "Number of indices: ",
    length(
      indices
    )
  )
  
  message(
    "Parallel workers: ",
    workers
  )
  
  message(
    "Depths per index: ",
    paste(
      vertical_positions,
      collapse = ", "
    )
  )
  
  parallel_results <- future.apply::future_lapply(
    indices,
    run_one_index,
    future.seed = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Combine site-depth summaries
  # ---------------------------------------------------------------------------
  
  summary_list <- lapply(
    parallel_results,
    
    function(
    result_i
    ) {
      
      getElement(
        result_i,
        "summary"
      )
    }
  )
  
  summary_list <- Filter(
    Negate(
      is.null
    ),
    summary_list
  )
  
  if (
    length(
      summary_list
    ) > 0L
  ) {
    
    combined_summary <- data.table::rbindlist(
      summary_list,
      use.names = TRUE,
      fill = TRUE
    )
    
    data.table::setorder(
      combined_summary,
      index,
      vertical_position
    )
    
  } else {
    
    combined_summary <-
      data.table::data.table()
  }
  
  # ---------------------------------------------------------------------------
  # Combine status
  # ---------------------------------------------------------------------------
  
  status_list <- lapply(
    parallel_results,
    
    function(
    result_i
    ) {
      
      getElement(
        result_i,
        "status"
      )
    }
  )
  
  combined_status <- data.table::rbindlist(
    status_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  data.table::setorder(
    combined_status,
    index,
    vertical_position
  )
  
  # ---------------------------------------------------------------------------
  # Save
  # ---------------------------------------------------------------------------
  
  dir.create(
    dirname(
      combined_output_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    dirname(
      status_output_file
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (
    nrow(
      combined_summary
    ) > 0L
  ) {
    
    data.table::fwrite(
      combined_summary,
      combined_output_file
    )
  }
  
  data.table::fwrite(
    combined_status,
    status_output_file
  )
  
  message(
    "\n========================================"
  )
  
  message(
    "PARALLEL RUN FINISHED"
  )
  
  message(
    "Indices attempted: ",
    length(
      indices
    )
  )
  
  message(
    "Successful index-depth runs: ",
    sum(
      combined_status$status ==
        "success"
    )
  )
  
  message(
    "Skipped because depth had no usable data: ",
    sum(
      combined_status$status ==
        "skipped_no_depth_data"
    )
  )
  
  message(
    "Other failures: ",
    sum(
      combined_status$status ==
        "failed"
    )
  )
  
  message(
    "========================================"
  )
  
  return(
    invisible(
      list(
        summary =
          combined_summary,
        
        status =
          combined_status
      )
    )
  )
}

# =============================================================================
# Permafrost shallow-soil temperature calibration
#
# Requires these generic SoilT helper functions to already exist:
#
#   causal_exponential_filter()
#   predict_daily_tsoil()
#   calc_metrics()
#   select_loyo_years()
#   read_sipnet_clim()
#   read_neon_shallow_soil_temperature()
#   prepare_tau_analysis_data()
#   NEON_NOMINAL_DEPTH_CM
#
# Permafrost model:
#
# T_eff =
#   a +
#   n_warm * max(Tair, 0) +
#   n_cold * min(Tair, 0)
#
# Tsoil(t) =
#   Tsoil(t-1) +
#   [1 - exp(-dt / tau)] *
#   [T_eff(t) - Tsoil(t-1)]
#
# Parameters:
#   a_C       = baseline soil-air thermal offset
#   n_warm    = warm-condition Tair -> Tsoil coupling
#   n_cold    = cold-condition Tair -> Tsoil coupling
#   tau_days  = thermal-memory timescale
# =============================================================================


# =============================================================================
# Dependency check
# =============================================================================

assert_permafrost_dependencies <- function() {
  
  required_functions <- c(
    "causal_exponential_filter",
    "predict_daily_tsoil",
    "calc_metrics",
    "select_loyo_years",
    "read_sipnet_clim",
    "read_neon_shallow_soil_temperature",
    "prepare_tau_analysis_data"
  )
  
  missing_functions <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
  ]
  
  if (
    length(missing_functions) > 0L
  ) {
    
    stop(
      "Missing required helper functions: ",
      paste(
        missing_functions,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

# =============================================================================
# Effective air-temperature forcing
# =============================================================================

permafrost_effective_tair <- function(
    tair_C,
    a_C,
    n_warm,
    n_cold
) {
  
  if (
    length(a_C) != 1L ||
    !is.finite(a_C)
  ) {
    
    stop(
      "`a_C` must be one finite value.",
      call. = FALSE
    )
  }
  
  if (
    length(n_warm) != 1L ||
    !is.finite(n_warm) ||
    n_warm < 0
  ) {
    
    stop(
      "`n_warm` must be one non-negative finite value.",
      call. = FALSE
    )
  }
  
  if (
    length(n_cold) != 1L ||
    !is.finite(n_cold) ||
    n_cold < 0
  ) {
    
    stop(
      "`n_cold` must be one non-negative finite value.",
      call. = FALSE
    )
  }
  
  effective_tair_C <-
    a_C +
    n_warm *
    pmax(
      tair_C,
      0
    ) +
    n_cold *
    pmin(
      tair_C,
      0
    )
  
  return(
    effective_tair_C
  )
}


# =============================================================================
# Subdaily process model
# =============================================================================

permafrost_tsoil_process_model <- function(
    tair_C,
    timestep_days,
    a_C,
    n_warm,
    n_cold,
    tau_days,
    initial_soil_temp_C = NULL
) {
  
  if (
    length(tau_days) != 1L ||
    !is.finite(tau_days) ||
    tau_days <= 0
  ) {
    
    stop(
      "`tau_days` must be one positive finite value.",
      call. = FALSE
    )
  }
  
  effective_tair_C <- permafrost_effective_tair(
    tair_C = tair_C,
    a_C = a_C,
    n_warm = n_warm,
    n_cold = n_cold
  )
  
  tsoil_C <- causal_exponential_filter(
    tair_C = effective_tair_C,
    timestep_days = timestep_days,
    tau_days = tau_days,
    initial_soil_temp_C = initial_soil_temp_C
  )
  
  return(
    tsoil_C
  )
}


# =============================================================================
# Daily prediction
# =============================================================================

predict_daily_permafrost_tsoil <- function(
    clim,
    a_C,
    n_warm,
    n_cold,
    tau_days
) {
  
  required_columns <- c(
    "date",
    "air_temp_C",
    "timestep_days"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(clim)
  )
  
  if (
    length(missing_columns) > 0L
  ) {
    
    stop(
      "Missing required climate columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  predicted_subdaily <- permafrost_tsoil_process_model(
    tair_C = clim$air_temp_C,
    timestep_days = clim$timestep_days,
    a_C = a_C,
    n_warm = n_warm,
    n_cold = n_cold,
    tau_days = tau_days
  )
  
  prediction <- data.table::data.table(
    date = clim$date,
    Tsoil_permafrost_C = predicted_subdaily
  )[
    ,
    .(
      Tsoil_permafrost_C =
        mean(
          Tsoil_permafrost_C,
          na.rm = TRUE
        )
    ),
    by = date
  ]
  
  data.table::setorder(
    prediction,
    date
  )
  
  return(
    prediction
  )
}


# =============================================================================
# Filtered warm/cold basis functions
# =============================================================================

build_permafrost_daily_basis <- function(
    clim,
    tau_days
) {
  
  required_columns <- c(
    "date",
    "air_temp_C",
    "timestep_days"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(clim)
  )
  
  if (
    length(missing_columns) > 0L
  ) {
    
    stop(
      "Missing required climate columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (
    length(tau_days) != 1L ||
    !is.finite(tau_days) ||
    tau_days <= 0
  ) {
    
    stop(
      "`tau_days` must be one positive finite value.",
      call. = FALSE
    )
  }
  
  warm_air_C <- pmax(
    clim$air_temp_C,
    0
  )
  
  cold_air_C <- pmin(
    clim$air_temp_C,
    0
  )
  
  warm_filtered <- causal_exponential_filter(
    tair_C = warm_air_C,
    timestep_days = clim$timestep_days,
    tau_days = tau_days
  )
  
  cold_filtered <- causal_exponential_filter(
    tair_C = cold_air_C,
    timestep_days = clim$timestep_days,
    tau_days = tau_days
  )
  
  basis <- data.table::data.table(
    date = clim$date,
    warm_basis = warm_filtered,
    cold_basis = cold_filtered
  )[
    ,
    .(
      warm_basis =
        mean(
          warm_basis,
          na.rm = TRUE
        ),
      
      cold_basis =
        mean(
          cold_basis,
          na.rm = TRUE
        )
    ),
    by = date
  ]
  
  data.table::setorder(
    basis,
    date
  )
  
  return(
    basis
  )
}


# =============================================================================
# Conditional likelihood at fixed tau
# =============================================================================

profile_permafrost_given_tau <- function(
    clim,
    obs_daily,
    tau_days,
    fit_dates = NULL,
    a_lower_C = -10,
    a_upper_C = 10,
    n_warm_lower = 0,
    n_warm_upper = 1.5,
    n_cold_lower = 0,
    n_cold_upper = 1.2,
    min_fit_days = 30L
) {
  
  required_obs_columns <- c(
    "date",
    "Tsoil_obs_C"
  )
  
  missing_obs_columns <- setdiff(
    required_obs_columns,
    names(obs_daily)
  )
  
  if (
    length(missing_obs_columns) > 0L
  ) {
    
    stop(
      "Missing observed SoilT columns: ",
      paste(
        missing_obs_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (
    !is.finite(a_lower_C) ||
    !is.finite(a_upper_C) ||
    a_upper_C <= a_lower_C
  ) {
    
    stop(
      "Invalid bounds for `a_C`.",
      call. = FALSE
    )
  }
  
  if (
    !is.finite(n_warm_lower) ||
    !is.finite(n_warm_upper) ||
    n_warm_lower < 0 ||
    n_warm_upper <= n_warm_lower
  ) {
    
    stop(
      "Invalid bounds for `n_warm`.",
      call. = FALSE
    )
  }
  
  if (
    !is.finite(n_cold_lower) ||
    !is.finite(n_cold_upper) ||
    n_cold_lower < 0 ||
    n_cold_upper <= n_cold_lower
  ) {
    
    stop(
      "Invalid bounds for `n_cold`.",
      call. = FALSE
    )
  }
  
  obs_use <- data.table::as.data.table(
    data.table::copy(
      obs_daily
    )
  )
  
  obs_use[
    ,
    date :=
      as.Date(
        date
      )
  ]
  
  if (
    !is.null(fit_dates)
  ) {
    
    fit_dates <- unique(
      as.Date(
        fit_dates
      )
    )
    
    fit_dates <- fit_dates[
      !is.na(
        fit_dates
      )
    ]
    
    obs_use <- obs_use[
      date %in%
        fit_dates
    ]
  }
  
  basis <- build_permafrost_daily_basis(
    clim = clim,
    tau_days = tau_days
  )
  
  dat <- merge(
    obs_use,
    basis,
    by = "date",
    all = FALSE
  )
  
  dat <- dat[
    is.finite(Tsoil_obs_C) &
      is.finite(warm_basis) &
      is.finite(cold_basis)
  ]
  
  n_fit <- nrow(
    dat
  )
  
  if (
    n_fit <
    as.integer(min_fit_days)
  ) {
    
    return(
      list(
        a_C = NA_real_,
        n_warm = NA_real_,
        n_cold = NA_real_,
        tau_days = tau_days,
        sigma_C = NA_real_,
        sse = Inf,
        nll = Inf,
        n = n_fit,
        coupling_convergence = NA_integer_
      )
    )
  }
  
  X <- cbind(
    intercept = 1,
    warm_basis = dat$warm_basis,
    cold_basis = dat$cold_basis
  )
  
  beta_start <- tryCatch(
    
    stats::lm.fit(
      x = X,
      y = dat$Tsoil_obs_C
    )$coefficients,
    
    error = function(e) {
      
      c(
        0,
        1,
        0.5
      )
    }
  )
  
  beta_start <- as.numeric(
    beta_start
  )
  
  if (
    length(beta_start) != 3L ||
    any(
      !is.finite(beta_start)
    )
  ) {
    
    beta_start <- c(
      0,
      1,
      0.5
    )
  }
  
  lower <- c(
    a_lower_C,
    n_warm_lower,
    n_cold_lower
  )
  
  upper <- c(
    a_upper_C,
    n_warm_upper,
    n_cold_upper
  )
  
  beta_start <- pmin(
    pmax(
      beta_start,
      lower
    ),
    upper
  )
  
  objective_sse <- function(
    par
  ) {
    
    predicted <-
      par[1L] +
      par[2L] *
      dat$warm_basis +
      par[3L] *
      dat$cold_basis
    
    residual <-
      dat$Tsoil_obs_C -
      predicted
    
    return(
      sum(
        residual^2
      )
    )
  }
  
  coupling_fit <- stats::optim(
    par = beta_start,
    fn = objective_sse,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      maxit = 500L
    )
  )
  
  beta_hat <- coupling_fit$par
  
  predicted <-
    beta_hat[1L] +
    beta_hat[2L] *
    dat$warm_basis +
    beta_hat[3L] *
    dat$cold_basis
  
  residual <-
    dat$Tsoil_obs_C -
    predicted
  
  sse <- sum(
    residual^2
  )
  
  sigma2_hat <-
    sse /
    n_fit
  
  if (
    !is.finite(sigma2_hat) ||
    sigma2_hat <= 0
  ) {
    
    return(
      list(
        a_C =
          as.numeric(
            beta_hat[1L]
          ),
        
        n_warm =
          as.numeric(
            beta_hat[2L]
          ),
        
        n_cold =
          as.numeric(
            beta_hat[3L]
          ),
        
        tau_days =
          as.numeric(
            tau_days
          ),
        
        sigma_C =
          NA_real_,
        
        sse =
          as.numeric(
            sse
          ),
        
        nll =
          Inf,
        
        n =
          n_fit,
        
        coupling_convergence =
          coupling_fit$convergence
      )
    )
  }
  
  sigma_hat <- sqrt(
    sigma2_hat
  )
  
  nll <-
    0.5 *
    n_fit *
    (
      log(
        2 *
          pi *
          sigma2_hat
      ) +
        1
    )
  
  return(
    list(
      a_C =
        as.numeric(
          beta_hat[1L]
        ),
      
      n_warm =
        as.numeric(
          beta_hat[2L]
        ),
      
      n_cold =
        as.numeric(
          beta_hat[3L]
        ),
      
      tau_days =
        as.numeric(
          tau_days
        ),
      
      sigma_C =
        as.numeric(
          sigma_hat
        ),
      
      sse =
        as.numeric(
          sse
        ),
      
      nll =
        as.numeric(
          nll
        ),
      
      n =
        n_fit,
      
      coupling_convergence =
        coupling_fit$convergence
    )
  )
}


# =============================================================================
# Full profiled Gaussian MLE
# =============================================================================

fit_permafrost_tsoil_mle <- function(
    clim,
    obs_daily,
    fit_dates = NULL,
    tau_lower_days = 0.5,
    tau_upper_days = 180,
    a_lower_C = -10,
    a_upper_C = 10,
    n_warm_lower = 0,
    n_warm_upper = 1.5,
    n_cold_lower = 0,
    n_cold_upper = 1.2,
    min_fit_days = 30L
) {
  
  if (
    length(tau_lower_days) != 1L ||
    length(tau_upper_days) != 1L ||
    !is.finite(tau_lower_days) ||
    !is.finite(tau_upper_days) ||
    tau_lower_days <= 0 ||
    tau_upper_days <= tau_lower_days
  ) {
    
    stop(
      paste0(
        "`tau_lower_days` and `tau_upper_days` must define ",
        "a positive increasing interval."
      ),
      call. = FALSE
    )
  }
  
  if (
    !is.null(fit_dates)
  ) {
    
    fit_dates <- unique(
      as.Date(
        fit_dates
      )
    )
    
    fit_dates <- fit_dates[
      !is.na(
        fit_dates
      )
    ]
    
    if (
      length(fit_dates) <
      as.integer(min_fit_days)
    ) {
      
      stop(
        "At least ",
        min_fit_days,
        " fitting dates are required.",
        call. = FALSE
      )
    }
  }
  
  objective_log_tau <- function(
    log_tau
  ) {
    
    tau_i <- exp(
      log_tau
    )
    
    fit_i <- profile_permafrost_given_tau(
      clim = clim,
      obs_daily = obs_daily,
      tau_days = tau_i,
      fit_dates = fit_dates,
      
      a_lower_C = a_lower_C,
      a_upper_C = a_upper_C,
      
      n_warm_lower = n_warm_lower,
      n_warm_upper = n_warm_upper,
      
      n_cold_lower = n_cold_lower,
      n_cold_upper = n_cold_upper,
      
      min_fit_days = min_fit_days
    )
    
    return(
      fit_i$nll
    )
  }
  
  tau_fit <- stats::optimize(
    f = objective_log_tau,
    
    interval = c(
      log(
        tau_lower_days
      ),
      log(
        tau_upper_days
      )
    ),
    
    tol = 1e-6
  )
  
  if (
    !is.finite(
      tau_fit$objective
    )
  ) {
    
    stop(
      "Permafrost SoilT MLE did not identify a finite likelihood.",
      call. = FALSE
    )
  }
  
  tau_hat <- exp(
    tau_fit$minimum
  )
  
  best <- profile_permafrost_given_tau(
    clim = clim,
    obs_daily = obs_daily,
    tau_days = tau_hat,
    fit_dates = fit_dates,
    
    a_lower_C = a_lower_C,
    a_upper_C = a_upper_C,
    
    n_warm_lower = n_warm_lower,
    n_warm_upper = n_warm_upper,
    
    n_cold_lower = n_cold_lower,
    n_cold_upper = n_cold_upper,
    
    min_fit_days = min_fit_days
  )
  
  tau_tol <- max(
    1e-6,
    1e-4 *
      (
        tau_upper_days -
          tau_lower_days
      )
  )
  
  a_tol <- max(
    1e-6,
    1e-4 *
      (
        a_upper_C -
          a_lower_C
      )
  )
  
  nw_tol <- max(
    1e-6,
    1e-4 *
      (
        n_warm_upper -
          n_warm_lower
      )
  )
  
  nc_tol <- max(
    1e-6,
    1e-4 *
      (
        n_cold_upper -
          n_cold_lower
      )
  )
  
  return(
    list(
      a_C =
        best$a_C,
      
      n_warm =
        best$n_warm,
      
      n_cold =
        best$n_cold,
      
      tau_days =
        best$tau_days,
      
      sigma_C =
        best$sigma_C,
      
      sse =
        best$sse,
      
      nll =
        best$nll,
      
      n_train =
        best$n,
      
      convergence =
        0L,
      
      coupling_convergence =
        best$coupling_convergence,
      
      tau_at_lower_bound =
        abs(
          best$tau_days -
            tau_lower_days
        ) <=
        tau_tol,
      
      tau_at_upper_bound =
        abs(
          best$tau_days -
            tau_upper_days
        ) <=
        tau_tol,
      
      a_at_bound =
        abs(
          best$a_C -
            a_lower_C
        ) <=
        a_tol ||
        abs(
          best$a_C -
            a_upper_C
        ) <=
        a_tol,
      
      n_warm_at_bound =
        abs(
          best$n_warm -
            n_warm_lower
        ) <=
        nw_tol ||
        abs(
          best$n_warm -
            n_warm_upper
        ) <=
        nw_tol,
      
      n_cold_at_bound =
        abs(
          best$n_cold -
            n_cold_lower
        ) <=
        nc_tol ||
        abs(
          best$n_cold -
            n_cold_upper
        ) <=
        nc_tol,
      
      optimizer =
        paste0(
          "Profiled Gaussian MLE: ",
          "bounded a_C/n_warm/n_cold conditional on tau; ",
          "stats::optimize on log(tau)"
        )
    )
  )
}


# =============================================================================
# Final corrected subdaily SoilT
# =============================================================================

make_permafrost_corrected_subdaily <- function(
    clim,
    final_fit
) {
  
  required_columns <- c(
    "timestamp",
    "date",
    "air_temp_C",
    "soil_temp_C",
    "timestep_days"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(clim)
  )
  
  if (
    length(missing_columns) > 0L
  ) {
    
    stop(
      "Missing climate columns: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  optimized_tsoil <- permafrost_tsoil_process_model(
    tair_C = clim$air_temp_C,
    timestep_days = clim$timestep_days,
    a_C = final_fit$a_C,
    n_warm = final_fit$n_warm,
    n_cold = final_fit$n_cold,
    tau_days = final_fit$tau_days
  )
  
  effective_tair <- permafrost_effective_tair(
    tair_C = clim$air_temp_C,
    a_C = final_fit$a_C,
    n_warm = final_fit$n_warm,
    n_cold = final_fit$n_cold
  )
  
  return(
    data.table::data.table(
      timestamp =
        clim$timestamp,
      
      date =
        clim$date,
      
      air_temp_C =
        clim$air_temp_C,
      
      effective_air_temp_C =
        effective_tair,
      
      soil_temp_current_C =
        clim$soil_temp_C,
      
      soil_temp_permafrost_opt_C =
        optimized_tsoil
    )
  )
}


# =============================================================================
# Leave-one-year-out validation
# =============================================================================

run_loyo_permafrost_tsoil_validation <- function(
    site,
    clim,
    obs_tsoil,
    analysis_data,
    tau_lower_days = 0.5,
    tau_upper_days = 180,
    a_lower_C = -10,
    a_upper_C = 10,
    n_warm_lower = 0,
    n_warm_upper = 1.5,
    n_cold_lower = 0,
    n_cold_upper = 1.2,
    min_days_per_year = 120L,
    min_days_per_season = 10L,
    min_train_days = 60L
) {
  
  year_info <- select_loyo_years(
    analysis_data = analysis_data,
    min_days_per_year = min_days_per_year,
    min_days_per_season = min_days_per_season
  )
  
  eligible_years <- year_info$eligible_years
  
  if (
    length(eligible_years) < 2L
  ) {
    
    stop(
      "LOYO requires at least two eligible years. Eligible years: ",
      paste(
        eligible_years,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # tau = 15 baseline
  # ---------------------------------------------------------------------------
  
  pred_tau15 <- predict_daily_tsoil(
    clim = clim,
    tau_days = 15,
    offset_C = 0
  )
  
  data.table::setnames(
    pred_tau15,
    "Tsoil_tau_C",
    "Tsoil_causal_tau15_C"
  )
  
  prediction_methods <- c(
    Current_SIPNET =
      "Tsoil_current_C",
    
    Causal_tau15 =
      "Tsoil_causal_tau15_C",
    
    Permafrost_MLE =
      "Tsoil_permafrost_opt_C"
  )
  
  parameter_rows <- list()
  metric_rows <- list()
  validation_rows <- list()
  
  # ---------------------------------------------------------------------------
  # LOYO
  # ---------------------------------------------------------------------------
  
  for (
    held_out_year_i in
    eligible_years
  ) {
    
    train_dates <- analysis_data[
      obs_year !=
        held_out_year_i,
      unique(
        date
      )
    ]
    
    if (
      length(train_dates) <
      as.integer(min_train_days)
    ) {
      
      warning(
        "Skipping held-out year ",
        held_out_year_i,
        ": only ",
        length(train_dates),
        " training days.",
        call. = FALSE
      )
      
      next
    }
    
    message(
      "LOYO ",
      site,
      " | held-out year = ",
      held_out_year_i,
      " | training days = ",
      length(train_dates)
    )
    
    fold_fit <- fit_permafrost_tsoil_mle(
      clim = clim,
      obs_daily = obs_tsoil,
      fit_dates = train_dates,
      
      tau_lower_days =
        tau_lower_days,
      
      tau_upper_days =
        tau_upper_days,
      
      a_lower_C =
        a_lower_C,
      
      a_upper_C =
        a_upper_C,
      
      n_warm_lower =
        n_warm_lower,
      
      n_warm_upper =
        n_warm_upper,
      
      n_cold_lower =
        n_cold_lower,
      
      n_cold_upper =
        n_cold_upper,
      
      min_fit_days =
        min_train_days
    )
    
    pred_optimized <- predict_daily_permafrost_tsoil(
      clim = clim,
      a_C = fold_fit$a_C,
      n_warm = fold_fit$n_warm,
      n_cold = fold_fit$n_cold,
      tau_days = fold_fit$tau_days
    )
    
    data.table::setnames(
      pred_optimized,
      "Tsoil_permafrost_C",
      "Tsoil_permafrost_opt_C"
    )
    
    fold_data <- Reduce(
      function(
    x,
    y
      ) {
        
        merge(
          x,
          y,
          by = "date",
          all = FALSE
        )
      },
    
    list(
      analysis_data,
      pred_tau15,
      pred_optimized
    )
    )
    
    test_data <- fold_data[
      obs_year ==
        held_out_year_i
    ]
    
    if (
      nrow(test_data) == 0L
    ) {
      
      next
    }
    
    test_data[
      ,
      held_out_year :=
        held_out_year_i
    ]
    
    validation_rows[
      length(validation_rows) +
        1L
    ] <- list(
      test_data
    )
    
    parameter_rows[
      length(parameter_rows) +
        1L
    ] <- list(
      data.table::data.table(
        held_out_year =
          held_out_year_i,
        
        a_C =
          fold_fit$a_C,
        
        n_warm =
          fold_fit$n_warm,
        
        n_cold =
          fold_fit$n_cold,
        
        tau_days =
          fold_fit$tau_days,
        
        sigma_C =
          fold_fit$sigma_C,
        
        nll =
          fold_fit$nll,
        
        n_train =
          fold_fit$n_train,
        
        tau_at_lower_bound =
          fold_fit$tau_at_lower_bound,
        
        tau_at_upper_bound =
          fold_fit$tau_at_upper_bound,
        
        a_at_bound =
          fold_fit$a_at_bound,
        
        n_warm_at_bound =
          fold_fit$n_warm_at_bound,
        
        n_cold_at_bound =
          fold_fit$n_cold_at_bound
      )
    )
    
    # -------------------------------------------------------------------------
    # Held-out metrics
    # -------------------------------------------------------------------------
    
    for (
      method_name in
      names(
        prediction_methods
      )
    ) {
      
      prediction_column <- unname(
        prediction_methods[
          method_name
        ]
      )
      
      predicted_values <- getElement(
        test_data,
        prediction_column
      )
      
      metrics <- calc_metrics(
        observed =
          test_data$Tsoil_obs_C,
        
        predicted =
          predicted_values
      )
      
      metrics[
        ,
        `:=`(
          held_out_year =
            held_out_year_i,
          
          method =
            method_name,
          
          fold_a_C =
            fold_fit$a_C,
          
          fold_n_warm =
            fold_fit$n_warm,
          
          fold_n_cold =
            fold_fit$n_cold,
          
          fold_tau_days =
            fold_fit$tau_days
        )
      ]
      
      metric_rows[
        length(metric_rows) +
          1L
      ] <- list(
        metrics
      )
    }
  }
  
  if (
    length(validation_rows) == 0L
  ) {
    
    stop(
      "No LOYO validation folds were completed.",
      call. = FALSE
    )
  }
  
  validation <- data.table::rbindlist(
    validation_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  parameter_folds <- data.table::rbindlist(
    parameter_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  fold_metrics <- data.table::rbindlist(
    metric_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Pooled out-of-sample metrics
  # ---------------------------------------------------------------------------
  
  pooled_rows <- list()
  
  for (
    method_name in
    names(
      prediction_methods
    )
  ) {
    
    prediction_column <- unname(
      prediction_methods[
        method_name
      ]
    )
    
    predicted_values <- getElement(
      validation,
      prediction_column
    )
    
    metrics <- calc_metrics(
      observed =
        validation$Tsoil_obs_C,
      
      predicted =
        predicted_values
    )
    
    metrics[
      ,
      `:=`(
        method =
          method_name,
        
        validation_scheme =
          "leave_one_well_observed_year_out",
        
        n_folds =
          data.table::uniqueN(
            validation$held_out_year
          )
      )
    ]
    
    pooled_rows[
      length(pooled_rows) +
        1L
    ] <- list(
      metrics
    )
  }
  
  pooled_metrics <- data.table::rbindlist(
    pooled_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  return(
    list(
      year_coverage =
        year_info$coverage,
      
      eligible_years =
        eligible_years,
      
      parameter_folds =
        parameter_folds,
      
      fold_metrics =
        fold_metrics,
      
      pooled_metrics =
        pooled_metrics,
      
      validation =
        validation,
      
      pred_tau15 =
        pred_tau15
    )
  )
}


# =============================================================================
# Complete workflow for one site and one depth
# =============================================================================

run_site_permafrost_tsoil_workflow <- function(
    config
) {
  
  assert_permafrost_dependencies()
  
  site <- config$site
  
  dir.create(
    config$output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  message(
    "Site: ",
    site,
    " | VER=",
    config$vertical_position,
    " | depth ~",
    config$nominal_depth_cm,
    " cm"
  )
  
  # ---------------------------------------------------------------------------
  # Read SIPNET forcing
  # ---------------------------------------------------------------------------
  
  clim <- read_sipnet_clim(
    config$clim_file
  )
  
  # ---------------------------------------------------------------------------
  # Read NEON SoilT
  # ---------------------------------------------------------------------------
  
  obs_tsoil <- read_neon_shallow_soil_temperature(
    site =
      site,
    
    env_dir =
      config$neon_env_dir,
    
    start_year =
      config$obs_start_year,
    
    end_year =
      config$obs_end_year,
    
    vertical_position =
      config$vertical_position,
    
    nominal_depth_cm =
      config$nominal_depth_cm,
    
    min_plots_per_timestamp =
      config$min_plots_per_timestamp,
    
    min_timestamps_per_day =
      config$min_tsoil_timestamps_per_day
  )
  
  # ---------------------------------------------------------------------------
  # Prepare daily overlap
  # ---------------------------------------------------------------------------
  
  analysis_data <- prepare_tau_analysis_data(
    clim =
      clim,
    
    obs_tsoil =
      obs_tsoil,
    
    obs_start_year =
      config$obs_start_year,
    
    obs_end_year =
      config$obs_end_year,
    
    warmup_days =
      config$warmup_days
  )
  
  # ---------------------------------------------------------------------------
  # LOYO validation
  # ---------------------------------------------------------------------------
  
  loyo <- run_loyo_permafrost_tsoil_validation(
    site =
      site,
    
    clim =
      clim,
    
    obs_tsoil =
      obs_tsoil,
    
    analysis_data =
      analysis_data,
    
    tau_lower_days =
      config$tau_lower_days,
    
    tau_upper_days =
      config$tau_upper_days,
    
    a_lower_C =
      config$a_lower_C,
    
    a_upper_C =
      config$a_upper_C,
    
    n_warm_lower =
      config$n_warm_lower,
    
    n_warm_upper =
      config$n_warm_upper,
    
    n_cold_lower =
      config$n_cold_lower,
    
    n_cold_upper =
      config$n_cold_upper,
    
    min_days_per_year =
      config$min_loyo_days_per_year,
    
    min_days_per_season =
      config$min_loyo_days_per_season,
    
    min_train_days =
      config$min_train_tsoil_days
  )
  
  # ---------------------------------------------------------------------------
  # Final all-data MLE
  #
  # No external "all_dates" object is needed.
  # ---------------------------------------------------------------------------
  
  final_fit <- fit_permafrost_tsoil_mle(
    clim =
      clim,
    
    obs_daily =
      obs_tsoil,
    
    fit_dates =
      analysis_data[
        ,
        unique(
          date
        )
      ],
    
    tau_lower_days =
      config$tau_lower_days,
    
    tau_upper_days =
      config$tau_upper_days,
    
    a_lower_C =
      config$a_lower_C,
    
    a_upper_C =
      config$a_upper_C,
    
    n_warm_lower =
      config$n_warm_lower,
    
    n_warm_upper =
      config$n_warm_upper,
    
    n_cold_lower =
      config$n_cold_lower,
    
    n_cold_upper =
      config$n_cold_upper,
    
    min_fit_days =
      config$min_train_tsoil_days
  )
  
  # ---------------------------------------------------------------------------
  # Final optimized subdaily SoilT
  # ---------------------------------------------------------------------------
  
  corrected_subdaily <- make_permafrost_corrected_subdaily(
    clim =
      clim,
    
    final_fit =
      final_fit
  )
  
  # ---------------------------------------------------------------------------
  # Pull pooled LOYO metrics
  # ---------------------------------------------------------------------------
  
  pooled_metrics <- loyo$pooled_metrics
  
  get_method_metric <- function(
    method_name,
    metric_name
  ) {
    
    method_rows <- pooled_metrics[
      method ==
        method_name
    ]
    
    if (
      nrow(method_rows) == 0L ||
      !metric_name %in%
      names(method_rows)
    ) {
      
      return(
        NA_real_
      )
    }
    
    metric_values <- getElement(
      method_rows,
      metric_name
    )
    
    return(
      as.numeric(
        metric_values[1L]
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Sensor depth metadata
  # ---------------------------------------------------------------------------
  
  zoffset_m <- attr(
    obs_tsoil,
    "zoffset_m_median"
  )
  
  observed_zoffset_m_median <- if (
    length(zoffset_m) == 1L &&
    is.finite(zoffset_m)
  ) {
    
    as.numeric(
      zoffset_m
    )
    
  } else {
    
    NA_real_
  }
  
  # ---------------------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------------------
  
  summary_table <- data.table::data.table(
    index =
      config$index,
    
    AmeriFlux_ID =
      config$AmeriFlux_ID,
    
    NEON_code =
      site,
    
    vertical_position =
      config$vertical_position,
    
    nominal_depth_cm =
      config$nominal_depth_cm,
    
    observed_zoffset_m_median =
      observed_zoffset_m_median,
    
    n_cv_years =
      nrow(
        loyo$parameter_folds
      ),
    
    # Final fitted parameters
    a_C =
      final_fit$a_C,
    
    n_warm =
      final_fit$n_warm,
    
    n_cold =
      final_fit$n_cold,
    
    tau_days =
      final_fit$tau_days,
    
    residual_sigma_C =
      final_fit$sigma_C,
    
    final_nll =
      final_fit$nll,
    
    # Boundary diagnostics
    tau_at_lower_bound =
      final_fit$tau_at_lower_bound,
    
    tau_at_upper_bound =
      final_fit$tau_at_upper_bound,
    
    a_at_bound =
      final_fit$a_at_bound,
    
    n_warm_at_bound =
      final_fit$n_warm_at_bound,
    
    n_cold_at_bound =
      final_fit$n_cold_at_bound,
    
    # Current PEcAn / SIPNET
    Current_SIPNET_LOYO_R2 =
      get_method_metric(
        "Current_SIPNET",
        "r2"
      ),
    
    Current_SIPNET_LOYO_RMSE_C =
      get_method_metric(
        "Current_SIPNET",
        "rmse"
      ),
    
    # Causal tau = 15
    Causal_tau15_LOYO_R2 =
      get_method_metric(
        "Causal_tau15",
        "r2"
      ),
    
    Causal_tau15_LOYO_RMSE_C =
      get_method_metric(
        "Causal_tau15",
        "rmse"
      ),
    
    # New permafrost model
    Permafrost_LOYO_R2 =
      get_method_metric(
        "Permafrost_MLE",
        "r2"
      ),
    
    Permafrost_LOYO_RMSE_C =
      get_method_metric(
        "Permafrost_MLE",
        "rmse"
      ),
    
    Permafrost_LOYO_MAE_C =
      get_method_metric(
        "Permafrost_MLE",
        "mae"
      ),
    
    Permafrost_LOYO_bias_C =
      get_method_metric(
        "Permafrost_MLE",
        "bias"
      ),
    
    Permafrost_LOYO_correlation =
      get_method_metric(
        "Permafrost_MLE",
        "correlation"
      )
  )
  
  # ---------------------------------------------------------------------------
  # Write outputs
  # ---------------------------------------------------------------------------
  
  data.table::fwrite(
    obs_tsoil,
    
    file.path(
      config$output_dir,
      paste0(
        site,
        "_NEON_tsoil_daily.csv"
      )
    )
  )
  
  data.table::fwrite(
    loyo$year_coverage,
    
    file.path(
      config$output_dir,
      paste0(
        site,
        "_LOYO_year_coverage.csv"
      )
    )
  )
  
  data.table::fwrite(
    loyo$validation,
    
    file.path(
      config$output_dir,
      paste0(
        site,
        "_permafrost_LOYO_predictions.csv"
      )
    )
  )
  
  data.table::fwrite(
    loyo$parameter_folds,
    
    file.path(
      config$output_dir,
      paste0(
        site,
        "_permafrost_LOYO_parameters.csv"
      )
    )
  )
  
  data.table::fwrite(
    loyo$fold_metrics,
    
    file.path(
      config$output_dir,
      paste0(
        site,
        "_permafrost_LOYO_fold_metrics.csv"
      )
    )
  )
  
  data.table::fwrite(
    loyo$pooled_metrics,
    
    file.path(
      config$output_dir,
      paste0(
        site,
        "_permafrost_LOYO_pooled_metrics.csv"
      )
    )
  )
  
  data.table::fwrite(
    corrected_subdaily,
    
    file.path(
      config$output_dir,
      paste0(
        site,
        "_permafrost_corrected_subdaily.csv"
      )
    )
  )
  
  summary_file <- file.path(
    config$output_dir,
    
    sprintf(
      "%s_index%d_VER%s_permafrost_summary.csv",
      
      site,
      
      as.integer(
        config$index
      ),
      
      as.character(
        config$vertical_position
      )
    )
  )
  
  data.table::fwrite(
    summary_table,
    summary_file
  )
  
  message(
    sprintf(
      paste0(
        "Final MLE | a=%.3f C | ",
        "n_warm=%.3f | ",
        "n_cold=%.3f | ",
        "tau=%.3f d | ",
        "LOYO R2=%.3f"
      ),
      
      final_fit$a_C,
      final_fit$n_warm,
      final_fit$n_cold,
      final_fit$tau_days,
      summary_table$Permafrost_LOYO_R2
    )
  )
  
  return(
    invisible(
      list(
        config =
          config,
        
        clim =
          clim,
        
        obs_tsoil =
          obs_tsoil,
        
        analysis_data =
          analysis_data,
        
        loyo =
          loyo,
        
        final_fit =
          final_fit,
        
        corrected_subdaily =
          corrected_subdaily,
        
        summary =
          summary_table,
        
        summary_file =
          summary_file
      )
    )
  )
}


# =============================================================================
# One lookup index and one depth
# =============================================================================

estimate_permafrost_tsoil_by_index <- function(
    index,
    lookup,
    vertical_position = "502",
    nominal_depth_cm = NULL,
    obs_start_year = 2017L,
    obs_end_year = 2024L,
    
    clim_root =
      "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024",
    
    clim_basename =
      "ERA5.1.2012-01-01.2024-12-31.clim",
    
    neon_env_dir =
      "/projectnb/dietzelab/jzobitz/02-NEON-sites/env-data",
    
    output_root =
      "/projectnb/dietzelab/guYANG/soilparam/permafrost_tsoil",
    
    # Parameter bounds
    tau_lower_days = 0.5,
    tau_upper_days = 180,
    
    a_lower_C = -10,
    a_upper_C = 10,
    
    n_warm_lower = 0,
    n_warm_upper = 1.5,
    
    n_cold_lower = 0,
    n_cold_upper = 1.2,
    
    # Data / LOYO
    min_plots_per_timestamp = 1L,
    min_tsoil_timestamps_per_day = 12L,
    
    warmup_days = 180L,
    
    min_loyo_days_per_year = 120L,
    min_loyo_days_per_season = 10L,
    
    min_train_tsoil_days = 60L
) {
  
  assert_permafrost_dependencies()
  
  target_index <- suppressWarnings(
    as.integer(
      index
    )
  )
  
  if (
    length(target_index) != 1L ||
    is.na(target_index)
  ) {
    
    stop(
      "`index` must be one valid integer.",
      call. = FALSE
    )
  }
  
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
  
  required_lookup_columns <- c(
    "index",
    "NEON_code"
  )
  
  missing_lookup_columns <- setdiff(
    required_lookup_columns,
    names(lookup_dt)
  )
  
  if (
    length(missing_lookup_columns) > 0L
  ) {
    
    stop(
      "lookup is missing: ",
      paste(
        missing_lookup_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  site_info <- lookup_dt[
    index ==
      target_index
  ]
  
  if (
    nrow(site_info) != 1L
  ) {
    
    stop(
      "Expected exactly one lookup row for index ",
      target_index,
      "; found ",
      nrow(site_info),
      ".",
      call. = FALSE
    )
  }
  
  neon_code <- as.character(
    site_info$NEON_code[1L]
  )
  
  if (
    is.na(neon_code) ||
    neon_code == ""
  ) {
    
    stop(
      "Index ",
      target_index,
      " has no valid NEON_code.",
      call. = FALSE
    )
  }
  
  ameriflux_id <- if (
    "AmeriFlux_ID" %in%
    names(site_info)
  ) {
    
    as.character(
      site_info$AmeriFlux_ID[1L]
    )
    
  } else {
    
    NA_character_
  }
  
  vertical_position <- as.character(
    vertical_position
  )[1L]
  
  if (
    is.null(nominal_depth_cm)
  ) {
    
    depth_map <- get(
      "NEON_NOMINAL_DEPTH_CM",
      inherits = TRUE
    )
    
    if (
      !vertical_position %in%
      names(depth_map)
    ) {
      
      stop(
        "No nominal depth is defined for vertical_position = ",
        vertical_position,
        ".",
        call. = FALSE
      )
    }
    
    nominal_depth_cm <- as.numeric(
      depth_map[
        vertical_position
      ]
    )[1L]
  }
  
  if (
    !is.finite(nominal_depth_cm)
  ) {
    
    stop(
      "`nominal_depth_cm` must be finite.",
      call. = FALSE
    )
  }
  
  clim_file <- file.path(
    clim_root,
    
    sprintf(
      "ERA5_%d_1",
      target_index
    ),
    
    clim_basename
  )
  
  if (
    !file.exists(clim_file)
  ) {
    
    stop(
      "CLIM file not found: ",
      clim_file,
      call. = FALSE
    )
  }
  
  output_dir <- file.path(
    output_root,
    
    sprintf(
      "%s_index%d_VER%s",
      neon_code,
      target_index,
      vertical_position
    )
  )
  
  config <- list(
    index =
      target_index,
    
    AmeriFlux_ID =
      ameriflux_id,
    
    site =
      neon_code,
    
    vertical_position =
      vertical_position,
    
    nominal_depth_cm =
      as.numeric(
        nominal_depth_cm
      ),
    
    obs_start_year =
      as.integer(
        obs_start_year
      ),
    
    obs_end_year =
      as.integer(
        obs_end_year
      ),
    
    clim_file =
      clim_file,
    
    neon_env_dir =
      neon_env_dir,
    
    output_dir =
      output_dir,
    
    tau_lower_days =
      as.numeric(
        tau_lower_days
      ),
    
    tau_upper_days =
      as.numeric(
        tau_upper_days
      ),
    
    a_lower_C =
      as.numeric(
        a_lower_C
      ),
    
    a_upper_C =
      as.numeric(
        a_upper_C
      ),
    
    n_warm_lower =
      as.numeric(
        n_warm_lower
      ),
    
    n_warm_upper =
      as.numeric(
        n_warm_upper
      ),
    
    n_cold_lower =
      as.numeric(
        n_cold_lower
      ),
    
    n_cold_upper =
      as.numeric(
        n_cold_upper
      ),
    
    min_plots_per_timestamp =
      as.integer(
        min_plots_per_timestamp
      ),
    
    min_tsoil_timestamps_per_day =
      as.integer(
        min_tsoil_timestamps_per_day
      ),
    
    warmup_days =
      as.integer(
        warmup_days
      ),
    
    min_loyo_days_per_year =
      as.integer(
        min_loyo_days_per_year
      ),
    
    min_loyo_days_per_season =
      as.integer(
        min_loyo_days_per_season
      ),
    
    min_train_tsoil_days =
      as.integer(
        min_train_tsoil_days
      )
  )
  
  results <- run_site_permafrost_tsoil_workflow(
    config =
      config
  )
  
  return(
    invisible(
      list(
        summary =
          results$summary,
        
        best_parameters =
          results$final_fit,
        
        results =
          results
      )
    )
  )
}


# =============================================================================
# Run all permafrost indices and depths in parallel
# =============================================================================

run_all_permafrost_indices_depths_parallel <- function(
    lookup,
    workers = 16L,
    
    vertical_positions = c(
      "501",
      "502",
      "503",
      "504"
    ),
    
    depth_map =
      NEON_NOMINAL_DEPTH_CM,
    
    output_root =
      "/projectnb/dietzelab/guYANG/soilparam/permafrost_tsoil",
    
    combined_output_file =
      file.path(
        output_root,
        "all_permafrost_all_depths_summary.csv"
      ),
    
    status_output_file =
      file.path(
        output_root,
        "all_permafrost_all_depths_status.csv"
      ),
    
    ...
) {
  
  assert_permafrost_dependencies()
  
  if (
    !requireNamespace(
      "future",
      quietly = TRUE
    )
  ) {
    
    stop(
      "Package 'future' is required.",
      call. = FALSE
    )
  }
  
  if (
    !requireNamespace(
      "future.apply",
      quietly = TRUE
    )
  ) {
    
    stop(
      "Package 'future.apply' is required.",
      call. = FALSE
    )
  }
  
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
  
  required_columns <- c(
    "index",
    "NEON_code"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(lookup_dt)
  )
  
  if (
    length(missing_columns) > 0L
  ) {
    
    stop(
      "lookup is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  vertical_positions <- as.character(
    vertical_positions
  )
  
  unknown_positions <- setdiff(
    vertical_positions,
    names(depth_map)
  )
  
  if (
    length(unknown_positions) > 0L
  ) {
    
    stop(
      "No nominal depth is defined for: ",
      paste(
        unknown_positions,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  indices <- sort(
    unique(
      as.integer(
        lookup_dt$index
      )
    )
  )
  
  indices <- indices[
    is.finite(
      indices
    )
  ]
  
  if (
    length(indices) == 0L
  ) {
    
    stop(
      "No valid indices found in lookup.",
      call. = FALSE
    )
  }
  
  dots <- list(
    ...
  )
  
  reserved_names <- intersect(
    names(dots),
    
    c(
      "index",
      "lookup",
      "vertical_position",
      "nominal_depth_cm",
      "output_root"
    )
  )
  
  if (
    length(reserved_names) > 0L
  ) {
    
    stop(
      "Do not pass these arguments through `...`: ",
      paste(
        reserved_names,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # One index per parallel worker
  # ---------------------------------------------------------------------------
  
  run_one_index <- function(
    index_i
  ) {
    
    result_rows <- list()
    status_rows <- list()
    
    site_rows <- lookup_dt[
      index ==
        index_i
    ]
    
    neon_code <- if (
      nrow(site_rows) == 1L
    ) {
      
      as.character(
        site_rows$NEON_code[1L]
      )
      
    } else {
      
      NA_character_
    }
    
    for (
      vertical_position_i in
      vertical_positions
    ) {
      
      depth_i_cm <- as.numeric(
        depth_map[
          vertical_position_i
        ]
      )[1L]
      
      message(
        "Index ",
        index_i,
        " | Site ",
        neon_code,
        " | VER=",
        vertical_position_i,
        " | depth ~",
        depth_i_cm,
        " cm"
      )
      
      fit_args <- c(
        list(
          index =
            index_i,
          
          lookup =
            lookup_dt,
          
          vertical_position =
            vertical_position_i,
          
          nominal_depth_cm =
            depth_i_cm,
          
          output_root =
            output_root
        ),
        
        dots
      )
      
      fit_i <- tryCatch(
        
        do.call(
          estimate_permafrost_tsoil_by_index,
          fit_args
        ),
        
        error = function(e) {
          
          e
        }
      )
      
      if (
        inherits(
          fit_i,
          "error"
        )
      ) {
        
        error_message <- conditionMessage(
          fit_i
        )
        
        no_depth_data <- grepl(
          paste0(
            "No usable NEON SoilT|",
            "No daily NEON SoilT|",
            "No usable SoilT|",
            "No NEON SoilT"
          ),
          
          error_message,
          
          ignore.case = TRUE
        )
        
        status_rows[
          length(status_rows) +
            1L
        ] <- list(
          data.table::data.table(
            index =
              index_i,
            
            NEON_code =
              neon_code,
            
            vertical_position =
              vertical_position_i,
            
            nominal_depth_cm =
              depth_i_cm,
            
            status =
              if (
                no_depth_data
              ) {
                
                "skipped_no_depth_data"
                
              } else {
                
                "failed"
              },
            
            message =
              error_message
          )
        )
        
        message(
          if (
            no_depth_data
          ) {
            
            "SKIP"
            
          } else {
            
            "FAILED"
          },
          
          " | index=",
          index_i,
          
          " | VER=",
          vertical_position_i,
          
          " | ",
          error_message
        )
        
        next
      }
      
      result_rows[
        length(result_rows) +
          1L
      ] <- list(
        fit_i$summary
      )
      
      status_rows[
        length(status_rows) +
          1L
      ] <- list(
        data.table::data.table(
          index =
            index_i,
          
          NEON_code =
            neon_code,
          
          vertical_position =
            vertical_position_i,
          
          nominal_depth_cm =
            depth_i_cm,
          
          status =
            "success",
          
          message =
            NA_character_
        )
      )
      
      message(
        sprintf(
          paste0(
            "SUCCESS | index=%d | VER=%s | ",
            "a=%.2f | ",
            "n_warm=%.3f | ",
            "n_cold=%.3f | ",
            "tau=%.2f d | ",
            "R2=%.3f"
          ),
          
          index_i,
          vertical_position_i,
          fit_i$best_parameters$a_C,
          fit_i$best_parameters$n_warm,
          fit_i$best_parameters$n_cold,
          fit_i$best_parameters$tau_days,
          fit_i$summary$Permafrost_LOYO_R2
        )
      )
    }
    
    index_summary <- if (
      length(result_rows) > 0L
    ) {
      
      data.table::rbindlist(
        result_rows,
        use.names = TRUE,
        fill = TRUE
      )
      
    } else {
      
      NULL
    }
    
    index_status <- if (
      length(status_rows) > 0L
    ) {
      
      data.table::rbindlist(
        status_rows,
        use.names = TRUE,
        fill = TRUE
      )
      
    } else {
      
      data.table::data.table(
        index =
          index_i,
        
        NEON_code =
          neon_code,
        
        vertical_position =
          NA_character_,
        
        nominal_depth_cm =
          NA_real_,
        
        status =
          "no_result",
        
        message =
          NA_character_
      )
    }
    
    return(
      list(
        summary =
          index_summary,
        
        status =
          index_status
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Parallel execution
  # ---------------------------------------------------------------------------
  
  future::plan(
    future::multisession,
    
    workers =
      as.integer(
        workers
      )
  )
  
  message(
    "Starting permafrost SoilT calibration..."
  )
  
  message(
    "Indices: ",
    length(indices)
  )
  
  message(
    "Workers: ",
    workers
  )
  
  message(
    "Depths: ",
    paste(
      vertical_positions,
      collapse = ", "
    )
  )
  
  parallel_results <- future.apply::future_lapply(
    indices,
    run_one_index,
    future.seed = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # Combine summaries
  # ---------------------------------------------------------------------------
  
  summary_list <- lapply(
    parallel_results,
    
    function(
    result_i
    ) {
      
      getElement(
        result_i,
        "summary"
      )
    }
  )
  
  summary_list <- Filter(
    Negate(
      is.null
    ),
    summary_list
  )
  
  combined_summary <- if (
    length(summary_list) > 0L
  ) {
    
    data.table::rbindlist(
      summary_list,
      use.names = TRUE,
      fill = TRUE
    )
    
  } else {
    
    data.table::data.table()
  }
  
  # ---------------------------------------------------------------------------
  # Combine status
  # ---------------------------------------------------------------------------
  
  status_list <- lapply(
    parallel_results,
    
    function(
    result_i
    ) {
      
      getElement(
        result_i,
        "status"
      )
    }
  )
  
  combined_status <- data.table::rbindlist(
    status_list,
    use.names = TRUE,
    fill = TRUE
  )
  
  if (
    nrow(combined_summary) > 0L
  ) {
    
    data.table::setorder(
      combined_summary,
      index,
      vertical_position
    )
  }
  
  data.table::setorder(
    combined_status,
    index,
    vertical_position
  )
  
  # ---------------------------------------------------------------------------
  # Write combined output
  # ---------------------------------------------------------------------------
  
  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  if (
    nrow(combined_summary) > 0L
  ) {
    
    data.table::fwrite(
      combined_summary,
      combined_output_file
    )
  }
  
  data.table::fwrite(
    combined_status,
    status_output_file
  )
  
  message(
    "========================================"
  )
  
  message(
    "PERMAFROST RUN FINISHED"
  )
  
  message(
    "Successful runs: ",
    sum(
      combined_status$status ==
        "success"
    )
  )
  
  message(
    "No-depth skips: ",
    sum(
      combined_status$status ==
        "skipped_no_depth_data"
    )
  )
  
  message(
    "Failures: ",
    sum(
      combined_status$status ==
        "failed"
    )
  )
  
  message(
    "========================================"
  )
  
  return(
    invisible(
      list(
        summary =
          combined_summary,
        
        status =
          combined_status
      )
    )
  )
}