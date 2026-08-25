# =============================================================================
# Standalone unified NEON soil-temperature MLE + LOYO workflow
# Version: 2026-08-25.3 (data.table logical-filter hotfix)
#
# Canonical workflow for both non-permafrost and permafrost NEON sites.
# Source only this file for SoilT MLE; it contains every required scientific
# function and avoids the legacy duplicate permafrost process-model definition.
#
# Core fixes included:
#   * perfect predictions return R2 = 1 / RMSE = 0 instead of invalid metrics;
#   * LOYO coverage is computed after ERA5 timestamp matching and warm-up;
#   * at least two valid held-out folds are required;
#   * fold-to-fold tau mean, SD, CV, and range quantify parameter stability;
#   * final all-year MLE is retained even when LOYO cannot be completed;
#   * failed joint optimizations are rejected;
#   * the caller's future plan is restored after parallel execution.
# =============================================================================
# =============================================================================
# Permafrost SoilT empirical MLE workflow
#
# Data inputs:
#   multi_site$soilT       : NEON 3-hour soil-temperature observations
#   era5_all$era5_mean     : ERA5 3-hour ensemble-mean meteorological forcing
#
# Model modes:
#   fit_tau_only = TRUE
#     a_C, n_warm, and n_cold are fixed; only tau_days is estimated by MLE.
#
#   fit_tau_only = FALSE
#     a_C, n_warm, n_cold, and tau_days are jointly estimated by bounded MLE.
#
# Top-level functions in this file:
#   1. permafrost_tsoil_process_model()
#   2. permafrost_tsoil_parameter_model()
#   3. permafrost_tsoil_data_model()
#   4. fit_permafrost_site_mle()
#   5. validate_permafrost_loyo()
#   6. fit_permafrost_all_sites()
#
# @author Yang Gu
# =============================================================================


#' Run the empirical permafrost soil-temperature process model
#'
#' Predicts shallow soil temperature from air temperature using a warm/cold
#' asymmetric effective thermal forcing followed by a first-order causal
#' thermal-memory filter.
#'
#' The effective forcing is
#'
#' \deqn{
#' T_{eff,t} = a + n_{warm}\max(T_{air,t},0)
#'             + n_{cold}\min(T_{air,t},0)
#' }
#'
#' and soil temperature evolves according to
#'
#' \deqn{
#' T_{soil,t} = T_{soil,t-1}
#'   + \left(1-\exp(-\Delta t/\tau)\right)
#'   \left(T_{eff,t}-T_{soil,t-1}\right).
#' }
#'
#' The function calculates the actual time interval between consecutive forcing
#' timestamps. Therefore, it can accommodate occasional forcing gaps without
#' assuming that every timestep is exactly three hours.
#'
#' @md
#' @param time POSIXct or POSIXct-coercible vector giving forcing timestamps.
#' @param tair_C Numeric vector of air temperature in degrees C.
#' @param a_C Baseline soil-air thermal offset in degrees C.
#' @param n_warm Non-negative warm-condition air-to-soil coupling coefficient.
#' @param n_cold Non-negative cold-condition air-to-soil coupling coefficient.
#' @param tau_days Positive thermal-memory time constant in days.
#' @param initial_soil_temp_C Optional initial soil temperature in degrees C.
#'   When `NULL`, the first effective air temperature is used.
#'
#' @return A data.table containing `time`, `Tair_C`, `Tair_effective_C`, and
#'   `Tsoil_pred_C`.
#'
#' @export
#' @author Yang Gu
permafrost_tsoil_process_model <- function(
    time,
    tair_C,
    a_C,
    n_warm,
    n_cold,
    tau_days,
    initial_soil_temp_C = NULL
) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
  
  if (length(time) != length(tair_C)) {
    stop("`time` and `tair_C` must have the same length.", call. = FALSE)
  }
  
  if (length(tair_C) == 0L) {
    stop("No air-temperature forcing was supplied.", call. = FALSE)
  }
  
  time <- as.POSIXct(
    time,
    tz = "UTC"
  )
  
  if (any(is.na(time))) {
    stop("`time` contains unparseable timestamps.", call. = FALSE)
  }
  
  tair_C <- as.numeric(
    tair_C
  )
  
  if (any(!is.finite(tair_C))) {
    stop("`tair_C` contains non-finite values.", call. = FALSE)
  }
  
  parameters <- c(
    a_C = as.numeric(a_C)[1L],
    n_warm = as.numeric(n_warm)[1L],
    n_cold = as.numeric(n_cold)[1L],
    tau_days = as.numeric(tau_days)[1L]
  )
  
  if (any(!is.finite(parameters))) {
    stop("All process-model parameters must be finite.", call. = FALSE)
  }
  
  if (
    parameters["n_warm"] < 0 ||
    parameters["n_cold"] < 0 ||
    parameters["tau_days"] <= 0
  ) {
    stop(
      paste0(
        "`n_warm` and `n_cold` must be non-negative, and ",
        "`tau_days` must be positive."
      ),
      call. = FALSE
    )
  }
  
  forcing <- data.table::data.table(
    time = time,
    Tair_C = tair_C
  )
  
  data.table::setorder(
    forcing,
    time
  )
  
  if (anyDuplicated(forcing$time) > 0L) {
    stop("Duplicate forcing timestamps are not allowed.", call. = FALSE)
  }
  
  effective_tair <-
    parameters["a_C"] +
    parameters["n_warm"] *
    pmax(
      forcing$Tair_C,
      0
    ) +
    parameters["n_cold"] *
    pmin(
      forcing$Tair_C,
      0
    )
  
  n_time <- nrow(
    forcing
  )
  
  timestep_days <- rep(
    NA_real_,
    n_time
  )
  
  if (n_time > 1L) {
    
    timestep_days[2:n_time] <-
      as.numeric(
        diff(forcing$time),
        units = "days"
      )
    
    valid_timestep <- timestep_days[
      is.finite(timestep_days) &
        timestep_days > 0
    ]
    
    if (length(valid_timestep) == 0L) {
      stop("Could not determine a positive forcing timestep.", call. = FALSE)
    }
    
    timestep_days[1L] <- stats::median(
      valid_timestep
    )
    
  } else {
    
    timestep_days[1L] <- 3 / 24
  }
  
  if (
    any(
      !is.finite(timestep_days) |
      timestep_days <= 0
    )
  ) {
    stop("Forcing timestamps must increase strictly through time.", call. = FALSE)
  }
  
  tsoil_pred <- rep(
    NA_real_,
    n_time
  )
  
  if (is.null(initial_soil_temp_C)) {
    
    tsoil_pred[1L] <- effective_tair[1L]
    
  } else {
    
    initial_value <- as.numeric(
      initial_soil_temp_C
    )[1L]
    
    if (!is.finite(initial_value)) {
      stop("`initial_soil_temp_C` must be finite when supplied.", call. = FALSE)
    }
    
    tsoil_pred[1L] <- initial_value
  }
  
  if (n_time > 1L) {
    
    for (time_i in 2:n_time) {
      
      alpha_i <-
        1 -
        exp(
          -timestep_days[time_i] /
            parameters["tau_days"]
        )
      
      tsoil_pred[time_i] <-
        tsoil_pred[time_i - 1L] +
        alpha_i *
        (
          effective_tair[time_i] -
            tsoil_pred[time_i - 1L]
        )
    }
  }
  
  return(
    data.table::data.table(
      time = forcing$time,
      Tair_C = forcing$Tair_C,
      Tair_effective_C = as.numeric(effective_tair),
      Tsoil_pred_C = as.numeric(tsoil_pred)
    )
  )
}


#' Define and validate permafrost SoilT model parameters
#'
#' Defines the parameter configuration used by the permafrost SoilT maximum-
#' likelihood workflow. The function supports two fitting modes.
#'
#' When `fit_tau_only = TRUE`, `a_C`, `n_warm`, and `n_cold` are fixed at the
#' supplied values and only `tau_days` is estimated. When
#' `fit_tau_only = FALSE`, all four process-model parameters are estimated
#' jointly within the supplied bounds.
#'
#' This is a parameter-constraint model, not a Bayesian prior model. No prior
#' probability density is imposed on any fitted parameter.
#'
#' @md
#' @param fit_tau_only Logical. If `TRUE`, only `tau_days` is estimated.
#' @param fixed_a_C Fixed baseline soil-air offset used in tau-only mode.
#' @param fixed_n_warm Fixed warm coupling coefficient used in tau-only mode.
#' @param fixed_n_cold Fixed cold coupling coefficient used in tau-only mode.
#' @param initial_parameters Optional named numeric vector with starting values
#'   for `a_C`, `n_warm`, `n_cold`, and `tau_days`. Used for four-parameter
#'   fitting. Its `tau_days` value is also used as the nominal starting tau in
#'   returned diagnostics for tau-only fitting.
#' @param a_bounds Numeric length-two lower and upper bounds for `a_C`.
#' @param n_warm_bounds Numeric length-two bounds for `n_warm`.
#' @param n_cold_bounds Numeric length-two bounds for `n_cold`.
#' @param tau_bounds Numeric length-two positive bounds for `tau_days`.
#'
#' @return A list containing fitting mode, starting parameters, fixed
#'   parameters, and optimizer bounds.
#'
#' @export
#' @author Yang Gu
permafrost_tsoil_parameter_model <- function(
    fit_tau_only = TRUE,
    fixed_a_C = 0,
    fixed_n_warm = 1,
    fixed_n_cold = 0.3,
    initial_parameters = NULL,
    a_bounds = c(-10, 10),
    n_warm_bounds = c(0, 1.5),
    n_cold_bounds = c(0, 1.2),
    tau_bounds = c(0.125, 180)
) {
  
  check_bounds <- function(
    bounds,
    name,
    positive_lower = FALSE
  ) {
    
    bounds <- as.numeric(
      bounds
    )
    
    if (
      length(bounds) != 2L ||
      any(!is.finite(bounds)) ||
      bounds[2L] <= bounds[1L]
    ) {
      stop(
        "`",
        name,
        "` must contain two finite increasing values.",
        call. = FALSE
      )
    }
    
    if (
      isTRUE(positive_lower) &&
      bounds[1L] <= 0
    ) {
      stop(
        "The lower `tau_bounds` value must be positive.",
        call. = FALSE
      )
    }
    
    return(
      bounds
    )
  }
  
  a_bounds <- check_bounds(
    a_bounds,
    "a_bounds"
  )
  
  n_warm_bounds <- check_bounds(
    n_warm_bounds,
    "n_warm_bounds"
  )
  
  n_cold_bounds <- check_bounds(
    n_cold_bounds,
    "n_cold_bounds"
  )
  
  tau_bounds <- check_bounds(
    tau_bounds,
    "tau_bounds",
    positive_lower = TRUE
  )
  
  if (
    n_warm_bounds[1L] < 0 ||
    n_cold_bounds[1L] < 0
  ) {
    stop(
      "Coupling-coefficient lower bounds cannot be negative.",
      call. = FALSE
    )
  }
  
  fixed_parameters <- c(
    a_C = as.numeric(fixed_a_C)[1L],
    n_warm = as.numeric(fixed_n_warm)[1L],
    n_cold = as.numeric(fixed_n_cold)[1L]
  )
  
  if (
    any(!is.finite(fixed_parameters)) ||
    fixed_parameters["n_warm"] < 0 ||
    fixed_parameters["n_cold"] < 0
  ) {
    stop(
      paste0(
        "Fixed parameters must be finite, and fixed coupling ",
        "coefficients must be non-negative."
      ),
      call. = FALSE
    )
  }
  
  if (is.null(initial_parameters)) {
    
    initial_parameters <- c(
      a_C = 0,
      n_warm = 1,
      n_cold = 0.5,
      tau_days = 15
    )
  }
  
  initial_names <- names(
    initial_parameters
  )
  
  initial_parameters <- as.numeric(
    initial_parameters
  )
  
  if (is.null(initial_names)) {
    
    if (length(initial_parameters) != 4L) {
      stop(
        "Unnamed `initial_parameters` must contain exactly four values.",
        call. = FALSE
      )
    }
    
    initial_names <- c(
      "a_C",
      "n_warm",
      "n_cold",
      "tau_days"
    )
  }
  
  names(initial_parameters) <- initial_names
  
  required_parameters <- c(
    "a_C",
    "n_warm",
    "n_cold",
    "tau_days"
  )
  
  if (
    !all(
      required_parameters %in%
      names(initial_parameters)
    )
  ) {
    stop(
      paste0(
        "`initial_parameters` must contain a_C, n_warm, ",
        "n_cold, and tau_days."
      ),
      call. = FALSE
    )
  }
  
  initial_parameters <- initial_parameters[
    required_parameters
  ]
  
  if (any(!is.finite(initial_parameters))) {
    stop("All initial parameter values must be finite.", call. = FALSE)
  }
  
  lower <- c(
    a_C = a_bounds[1L],
    n_warm = n_warm_bounds[1L],
    n_cold = n_cold_bounds[1L],
    tau_days = tau_bounds[1L]
  )
  
  upper <- c(
    a_C = a_bounds[2L],
    n_warm = n_warm_bounds[2L],
    n_cold = n_cold_bounds[2L],
    tau_days = tau_bounds[2L]
  )
  
  initial_parameters <- pmin(
    pmax(
      initial_parameters,
      lower
    ),
    upper
  )
  
  if (isTRUE(fit_tau_only)) {
    
    initial_parameters["a_C"] <- fixed_parameters["a_C"]
    initial_parameters["n_warm"] <- fixed_parameters["n_warm"]
    initial_parameters["n_cold"] <- fixed_parameters["n_cold"]
  }
  
  return(
    list(
      fit_tau_only = isTRUE(fit_tau_only),
      fit_mode = if (isTRUE(fit_tau_only)) {
        "tau_only"
      } else {
        "four_parameter"
      },
      parameters = initial_parameters,
      fixed_parameters = fixed_parameters,
      lower = lower,
      upper = upper,
      tau_bounds = tau_bounds
    )
  )
}


#' Evaluate the Gaussian permafrost SoilT observation data model
#'
#' Evaluates a Gaussian observation model for paired NEON and predicted soil
#' temperatures. Residual standard deviation is profiled analytically at each
#' candidate process-model parameter set.
#'
#' The observation model is
#'
#' \deqn{
#' T_{soil,obs} = T_{soil,pred} + \epsilon,
#' \qquad \epsilon \sim Normal(0,\sigma^2),
#' }
#'
#' with
#'
#' \deqn{
#' \hat{\sigma}^2 = SSE/n.
#' }
#'
#' @md
#' @param observed_C Numeric vector of observed NEON soil temperatures.
#' @param predicted_C Numeric vector of process-model predictions.
#' @param min_observations Minimum number of finite paired observations required
#'   to evaluate the likelihood.
#'
#' @return A list containing `n`, `sigma_C`, `sse`, `nll`, residuals,
#'   predictive R2, RMSE, MAE, bias, and Pearson correlation.
#'
#' @export
#' @author Yang Gu
permafrost_tsoil_data_model <- function(
    observed_C,
    predicted_C,
    min_observations = 30L
) {
  min_observations <- suppressWarnings(
    as.integer(min_observations)[1L]
  )

  if (is.na(min_observations) || min_observations < 1L) {
    stop("`min_observations` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(observed_C) || !is.numeric(predicted_C)) {
    stop("`observed_C` and `predicted_C` must be numeric.", call. = FALSE)
  }
  
  if (length(observed_C) != length(predicted_C)) {
    stop(
      "`observed_C` and `predicted_C` must have the same length.",
      call. = FALSE
    )
  }
  
  valid <-
    is.finite(observed_C) &
    is.finite(predicted_C)
  
  observed <- as.numeric(
    observed_C[valid]
  )
  
  predicted <- as.numeric(
    predicted_C[valid]
  )
  
  n_obs <- length(
    observed
  )
  
  if (n_obs < min_observations) {
    
    return(
      list(
        n = n_obs,
        sigma_C = NA_real_,
        sse = Inf,
        nll = Inf,
        residual = numeric(),
        r2 = NA_real_,
        rmse_C = NA_real_,
        mae_C = NA_real_,
        bias_C = NA_real_,
        correlation = NA_real_
      )
    )
  }
  
  residual <- observed - predicted
  
  sse <- sum(
    residual^2
  )
  
  sigma2_mle <-
    sse /
    n_obs
  
  if (
    !is.finite(sigma2_mle) ||
    sigma2_mle < 0
  ) {
    
    return(
      list(
        n = n_obs,
        sigma_C = NA_real_,
        sse = sse,
        nll = Inf,
        residual = residual,
        r2 = NA_real_,
        rmse_C = NA_real_,
        mae_C = NA_real_,
        bias_C = NA_real_,
        correlation = NA_real_
      )
    )
  }
  
  # A zero SSE is a valid perfect prediction. Use a numerical variance floor
  # only for the Gaussian likelihood while retaining sigma_C = 0 and exact
  # predictive metrics.
  sigma_C <- sqrt(
    sigma2_mle
  )
  
  sigma2_for_likelihood <- max(
    sigma2_mle,
    .Machine$double.eps
  )
  
  nll <-
    0.5 *
    (
      n_obs *
        log(
          2 *
            pi *
            sigma2_for_likelihood
        ) +
        sse /
          sigma2_for_likelihood
    )
  
  denominator <- sum(
    (
      observed -
        mean(observed)
    )^2
  )
  
  r2 <- if (
    is.finite(denominator) &&
    denominator > 0
  ) {
    
    1 -
      sse /
      denominator
    
  } else {
    
    NA_real_
  }
  
  sd_observed <- if (n_obs > 1L) {
    stats::sd(observed)
  } else {
    NA_real_
  }
  
  sd_predicted <- if (n_obs > 1L) {
    stats::sd(predicted)
  } else {
    NA_real_
  }
  
  correlation <- if (
    is.finite(sd_observed) &&
    is.finite(sd_predicted) &&
    sd_observed > 0 &&
    sd_predicted > 0
  ) {
    
    stats::cor(
      observed,
      predicted
    )
    
  } else {
    
    NA_real_
  }
  
  return(
    list(
      n = n_obs,
      sigma_C = sigma_C,
      sse = sse,
      nll = nll,
      residual = residual,
      r2 = r2,
      rmse_C = sqrt(mean(residual^2)),
      mae_C = mean(abs(residual)),
      bias_C = mean(predicted - observed),
      correlation = correlation
    )
  )
}


#' Fit the permafrost SoilT model for one site and one NEON depth
#'
#' Matches one model-site index to pre-extracted NEON soil temperature and ERA5
#' ensemble-mean air temperature, runs the empirical permafrost process model,
#' and estimates process-model parameters by Gaussian maximum likelihood.
#'
#' With `fit_tau_only = TRUE`, `a_C`, `n_warm`, and `n_cold` are fixed and only
#' `tau_days` is estimated. The optimization is one-dimensional in `log(tau)`.
#' With `fit_tau_only = FALSE`, all four process parameters are jointly
#' estimated using bounded L-BFGS-B optimization.
#'
#' The process model always runs continuously across the complete requested
#' forcing period. `fit_years` restricts only which NEON observations contribute
#' to the likelihood. This prevents artificial process-state resets during
#' leave-one-year-out validation.
#'
#' @md
#' @param lookup Site lookup table containing `index`; `NEON_code` and
#'   `AmeriFlux_ID` are retained when available.
#' @param index Integer model-site index to fit.
#' @param start_year First calendar year included in the forcing period.
#' @param end_year Last calendar year included in the forcing period.
#' @param multi_site Object returned by `get_soil_neon_data_multi_site()`.
#' @param era5_all Object returned by `get_soil_era5_data_multi_site()`.
#' @param vertical_position Character NEON soil-temperature vertical position.
#' @param fit_years Optional integer vector of observation years used in the
#'   likelihood. `NULL` uses every available observation year.
#' @param fit_tau_only Logical. If `TRUE`, only `tau_days` is estimated.
#' @param fixed_a_C Fixed `a_C` used when `fit_tau_only = TRUE`.
#' @param fixed_n_warm Fixed `n_warm` used when `fit_tau_only = TRUE`.
#' @param fixed_n_cold Fixed `n_cold` used when `fit_tau_only = TRUE`.
#' @param warmup_days Number of days at the beginning of the forcing period
#'   excluded from the likelihood.
#' @param min_observations Minimum number of matched 3-hour observations needed
#'   for fitting.
#' @param a_bounds Bounds for `a_C` in four-parameter mode.
#' @param n_warm_bounds Bounds for `n_warm` in four-parameter mode.
#' @param n_cold_bounds Bounds for `n_cold` in four-parameter mode.
#' @param tau_bounds Positive bounds for `tau_days`.
#' @param initial_parameters Optional named starting values for `a_C`,
#'   `n_warm`, `n_cold`, and `tau_days`.
#' @param maxit Maximum L-BFGS-B iterations in four-parameter mode.
#'
#' @return A list containing site metadata, fitting mode, fitted process
#'   parameters, likelihood diagnostics, matched observations, and the complete
#'   forcing-period prediction.
#'
#' @export
#' @author Yang Gu
fit_permafrost_site_mle <- function(
    lookup,
    index,
    start_year,
    end_year,
    multi_site,
    era5_all,
    vertical_position = "501",
    fit_years = NULL,
    fit_tau_only = TRUE,
    fixed_a_C = 0,
    fixed_n_warm = 1,
    fixed_n_cold = 0.3,
    warmup_days = 180L,
    min_observations = 100L,
    a_bounds = c(-10, 10),
    n_warm_bounds = c(0, 1.5),
    n_cold_bounds = c(0, 1.2),
    tau_bounds = c(0.125, 180),
    initial_parameters = NULL,
    maxit = 1000L
) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
  
  site_index <- as.integer(
    index
  )[1L]
  
  start_year <- as.integer(
    start_year
  )[1L]
  
  end_year <- as.integer(
    end_year
  )[1L]
  
  vertical_position <- as.character(
    vertical_position
  )[1L]

  warmup_days <- suppressWarnings(
    as.numeric(warmup_days)[1L]
  )

  min_observations <- suppressWarnings(
    as.integer(min_observations)[1L]
  )

  maxit <- suppressWarnings(
    as.integer(maxit)[1L]
  )
  
  if (!is.finite(site_index)) {
    stop("`index` must be one finite integer.", call. = FALSE)
  }
  
  if (
    !is.finite(start_year) ||
    !is.finite(end_year) ||
    end_year < start_year
  ) {
    stop("Invalid `start_year` / `end_year`.", call. = FALSE)
  }

  if (is.na(vertical_position) || !nzchar(vertical_position)) {
    stop("`vertical_position` must be one non-empty value.", call. = FALSE)
  }

  if (!is.finite(warmup_days) || warmup_days < 0) {
    stop("`warmup_days` must be one non-negative value.", call. = FALSE)
  }

  if (is.na(min_observations) || min_observations < 1L) {
    stop("`min_observations` must be a positive integer.", call. = FALSE)
  }

  if (is.na(maxit) || maxit < 1L) {
    stop("`maxit` must be a positive integer.", call. = FALSE)
  }
  
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
  
  if (!"index" %in% names(lookup_dt)) {
    stop("`lookup` must contain `index`.", call. = FALSE)
  }
  
  site_row <- lookup_dt[
    index == site_index
  ]
  
  if (nrow(site_row) == 0L) {
    stop(
      "Index ",
      site_index,
      " was not found in lookup.",
      call. = FALSE
    )
  }
  
  site_row <- site_row[1L]
  
  soilT_all <- getElement(
    multi_site,
    "soilT"
  )
  
  era5_mean_all <- getElement(
    era5_all,
    "era5_mean"
  )
  
  if (is.null(soilT_all)) {
    stop("`multi_site` does not contain `soilT`.", call. = FALSE)
  }
  
  if (is.null(era5_mean_all)) {
    stop("`era5_all` does not contain `era5_mean`.", call. = FALSE)
  }
  
  soilT_all <- data.table::as.data.table(
    data.table::copy(
      soilT_all
    )
  )
  
  era5_mean_all <- data.table::as.data.table(
    data.table::copy(
      era5_mean_all
    )
  )
  
  required_soil_columns <- c(
    "time",
    "index",
    "verticalPosition",
    "SoilT"
  )
  
  missing_soil_columns <- setdiff(
    required_soil_columns,
    names(soilT_all)
  )
  
  if (length(missing_soil_columns) > 0L) {
    stop(
      "`multi_site$soilT` is missing: ",
      paste(
        missing_soil_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  if (
    !all(
      c(
        "time",
        "index"
      ) %in%
      names(era5_mean_all)
    )
  ) {
    stop(
      "`era5_all$era5_mean` must contain `time` and `index`.",
      call. = FALSE
    )
  }
  
  if (!"AirT_C" %in% names(era5_mean_all)) {
    
    if ("air_temperature" %in% names(era5_mean_all)) {
      
      era5_mean_all[
        ,
        AirT_C :=
          as.numeric(air_temperature) -
          273.15
      ]
      
    } else {
      
      stop(
        paste0(
          "`era5_all$era5_mean` must contain `AirT_C` ",
          "or `air_temperature`."
        ),
        call. = FALSE
      )
    }
  }
  
  soilT_all[
    ,
    time :=
      as.POSIXct(
        time,
        tz = "UTC"
      )
  ]
  
  era5_mean_all[
    ,
    time :=
      as.POSIXct(
        time,
        tz = "UTC"
      )
  ]
  
  start_time <- as.POSIXct(
    sprintf(
      "%04d-01-01 00:00:00",
      start_year
    ),
    tz = "UTC"
  )
  
  end_time <- as.POSIXct(
    sprintf(
      "%04d-01-01 00:00:00",
      end_year + 1L
    ),
    tz = "UTC"
  )
  
  soil_site <- soilT_all[
    index == site_index &
      as.character(verticalPosition) == vertical_position &
      time >= start_time &
      time < end_time &
      is.finite(SoilT)
  ]
  
  if (nrow(soil_site) == 0L) {
    
    available_depths <- unique(
      as.character(
        soilT_all[
          index == site_index,
          verticalPosition
        ]
      )
    )
    
    available_depths <- available_depths[
      !is.na(available_depths) &
        available_depths != ""
    ]
    
    stop(
      "No SoilT observations for index ",
      site_index,
      " at verticalPosition=",
      vertical_position,
      ". Available depths: ",
      paste(
        sort(available_depths),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  has_zoffset <-
    "mean_zOffset" %in%
    names(soil_site)
  
  if (has_zoffset) {
    
    soil_site <- soil_site[
      ,
      .(
        SoilT_obs_C =
          mean(
            SoilT,
            na.rm = TRUE
          ),
        mean_zOffset = {
          z_use <- mean_zOffset[
            is.finite(mean_zOffset)
          ]
          
          if (length(z_use) > 0L) {
            mean(z_use)
          } else {
            NA_real_
          }
        }
      ),
      by = time
    ]
    
  } else {
    
    soil_site <- soil_site[
      ,
      .(
        SoilT_obs_C =
          mean(
            SoilT,
            na.rm = TRUE
          ),
        mean_zOffset =
          NA_real_
      ),
      by = time
    ]
  }
  
  met_site <- era5_mean_all[
    index == site_index &
      time >= start_time &
      time < end_time &
      is.finite(AirT_C),
    .(
      AirT_C =
        mean(
          AirT_C,
          na.rm = TRUE
        )
    ),
    by = time
  ]
  
  data.table::setorder(
    met_site,
    time
  )
  
  if (nrow(met_site) < 2L) {
    stop(
      "Insufficient ERA5 air-temperature forcing for index ",
      site_index,
      ".",
      call. = FALSE
    )
  }
  
  warmup_end <-
    min(met_site$time) +
    as.numeric(warmup_days) *
    86400
  
  soil_fit <- soil_site[
    time >= warmup_end
  ]
  
  soil_fit[
    ,
    obs_year :=
      as.integer(
        format(
          time,
          "%Y",
          tz = "UTC"
        )
      )
  ]
  
  if (!is.null(fit_years)) {
    
    fit_years <- unique(
      as.integer(
        fit_years
      )
    )
    
    fit_years <- fit_years[
      is.finite(fit_years)
    ]
    
    soil_fit <- soil_fit[
      obs_year %in%
        fit_years
    ]
  }
  
  forcing_row <- match(
    soil_fit$time,
    met_site$time
  )
  
  soil_fit[
    ,
    forcing_row :=
      forcing_row
  ]
  
  soil_fit <- soil_fit[
    !is.na(forcing_row)
  ]
  
  if (
    nrow(soil_fit) <
    as.integer(min_observations)
  ) {
    stop(
      "Only ",
      nrow(soil_fit),
      " matched SoilT observations are available for fitting; minimum is ",
      min_observations,
      ".",
      call. = FALSE
    )
  }
  
  parameter_spec <- permafrost_tsoil_parameter_model(
    fit_tau_only = fit_tau_only,
    fixed_a_C = fixed_a_C,
    fixed_n_warm = fixed_n_warm,
    fixed_n_cold = fixed_n_cold,
    initial_parameters = initial_parameters,
    a_bounds = a_bounds,
    n_warm_bounds = n_warm_bounds,
    n_cold_bounds = n_cold_bounds,
    tau_bounds = tau_bounds
  )
  
  fit_tau_only <- getElement(
    parameter_spec,
    "fit_tau_only"
  )
  
  fit_mode <- getElement(
    parameter_spec,
    "fit_mode"
  )
  
  par_start <- getElement(
    parameter_spec,
    "parameters"
  )
  
  lower <- getElement(
    parameter_spec,
    "lower"
  )
  
  upper <- getElement(
    parameter_spec,
    "upper"
  )
  
  fixed_parameters <- getElement(
    parameter_spec,
    "fixed_parameters"
  )

  tau_bounds <- getElement(
    parameter_spec,
    "tau_bounds"
  )
  
  evaluate_process_parameters <- function(
    a_i,
    n_warm_i,
    n_cold_i,
    tau_i
  ) {
    
    process_prediction <- permafrost_tsoil_process_model(
      time = met_site$time,
      tair_C = met_site$AirT_C,
      a_C = a_i,
      n_warm = n_warm_i,
      n_cold = n_cold_i,
      tau_days = tau_i
    )
    
    predicted_full <- getElement(
      process_prediction,
      "Tsoil_pred_C"
    )
    
    predicted_fit <- predicted_full[
      soil_fit$forcing_row
    ]
    
    likelihood <- permafrost_tsoil_data_model(
      observed_C = soil_fit$SoilT_obs_C,
      predicted_C = predicted_fit,
      min_observations = min_observations
    )
    
    return(
      list(
        prediction = process_prediction,
        predicted_fit = predicted_fit,
        likelihood = likelihood
      )
    )
  }
  
  if (isTRUE(fit_tau_only)) {
    
    objective_log_tau <- function(
    log_tau
    ) {
      
      tau_i <- exp(
        log_tau
      )
      
      evaluation <- evaluate_process_parameters(
        a_i = fixed_parameters["a_C"],
        n_warm_i = fixed_parameters["n_warm"],
        n_cold_i = fixed_parameters["n_cold"],
        tau_i = tau_i
      )
      
      likelihood <- getElement(
        evaluation,
        "likelihood"
      )
      
      return(
        getElement(
          likelihood,
          "nll"
        )
      )
    }
    
    tau_fit <- stats::optimize(
      f = objective_log_tau,
      interval = log(tau_bounds),
      tol = 1e-7
    )
    
    if (!is.finite(tau_fit$objective)) {
      stop(
        "Tau-only optimization did not identify a finite likelihood.",
        call. = FALSE
      )
    }
    
    tau_hat <- exp(
      tau_fit$minimum
    )
    
    # IMPORTANT:
    # Drop the names on scalar values before rebuilding the named vector.
    # Otherwise R produces names such as `a_C.a_C`, causing later lookup of
    # `best_par["a_C"]` to return NA.
    best_par <- c(
      a_C = as.numeric(fixed_parameters["a_C"]),
      n_warm = as.numeric(fixed_parameters["n_warm"]),
      n_cold = as.numeric(fixed_parameters["n_cold"]),
      tau_days = as.numeric(tau_hat)
    )
    
    optimizer_convergence <- 0L
    optimizer_message <-
      "stats::optimize on log(tau); a_C, n_warm, and n_cold fixed"
    
  } else {
    
    objective_joint <- function(
    par
    ) {
      
      evaluation <- evaluate_process_parameters(
        a_i = par[1L],
        n_warm_i = par[2L],
        n_cold_i = par[3L],
        tau_i = par[4L]
      )
      
      likelihood <- getElement(
        evaluation,
        "likelihood"
      )
      
      return(
        getElement(
          likelihood,
          "nll"
        )
      )
    }
    
    fit <- stats::optim(
      par = par_start,
      fn = objective_joint,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(
        maxit = as.integer(maxit),
        parscale = c(
          5,
          0.5,
          0.5,
          20
        )
      )
    )
    
    if (
      !is.finite(fit$value) ||
      !identical(as.integer(fit$convergence), 0L)
    ) {
      stop(
        "Joint L-BFGS-B optimization failed to converge: ",
        if (is.null(fit$message)) {
          paste0("convergence code ", fit$convergence)
        } else {
          as.character(fit$message)
        },
        call. = FALSE
      )
    }
    
    best_par <- fit$par
    
    names(best_par) <- c(
      "a_C",
      "n_warm",
      "n_cold",
      "tau_days"
    )
    
    optimizer_convergence <- fit$convergence
    
    optimizer_message <- if (
      is.null(fit$message)
    ) {
      NA_character_
    } else {
      as.character(fit$message)
    }
  }
  
  final_evaluation <- evaluate_process_parameters(
    a_i = best_par["a_C"],
    n_warm_i = best_par["n_warm"],
    n_cold_i = best_par["n_cold"],
    tau_i = best_par["tau_days"]
  )
  
  final_prediction <- getElement(
    final_evaluation,
    "prediction"
  )
  
  final_fit_prediction <- getElement(
    final_evaluation,
    "predicted_fit"
  )
  
  final_data_model <- getElement(
    final_evaluation,
    "likelihood"
  )
  
  matched_fit <- data.table::copy(
    soil_fit
  )
  
  matched_fit[
    ,
    Tsoil_pred_C :=
      final_fit_prediction
  ]
  
  neon_code <- if (
    "NEON_code" %in%
    names(site_row)
  ) {
    as.character(
      site_row$NEON_code[1L]
    )
  } else {
    NA_character_
  }
  
  ameriflux_id <- if (
    "AmeriFlux_ID" %in%
    names(site_row)
  ) {
    as.character(
      site_row$AmeriFlux_ID[1L]
    )
  } else {
    NA_character_
  }
  
  mean_zoffset <- if (
    any(
      is.finite(
        matched_fit$mean_zOffset
      )
    )
  ) {
    
    mean(
      matched_fit$mean_zOffset[
        is.finite(matched_fit$mean_zOffset)
      ],
      na.rm = TRUE
    )
    
  } else {
    
    NA_real_
  }
  
  tau_tolerance <- max(
    1e-6,
    1e-4 *
      (
        tau_bounds[2L] -
          tau_bounds[1L]
      )
  )
  
  return(
    list(
      index = site_index,
      NEON_code = neon_code,
      AmeriFlux_ID = ameriflux_id,
      verticalPosition = vertical_position,
      mean_zOffset = mean_zoffset,
      start_year = start_year,
      end_year = end_year,
      fit_years = fit_years,
      fit_tau_only = isTRUE(fit_tau_only),
      fit_mode = fit_mode,
      fixed_a_C = as.numeric(fixed_parameters["a_C"]),
      fixed_n_warm = as.numeric(fixed_parameters["n_warm"]),
      fixed_n_cold = as.numeric(fixed_parameters["n_cold"]),
      parameters = best_par,
      sigma_C = getElement(final_data_model, "sigma_C"),
      nll = getElement(final_data_model, "nll"),
      n_fit = getElement(final_data_model, "n"),
      training_r2 = getElement(final_data_model, "r2"),
      training_rmse_C = getElement(final_data_model, "rmse_C"),
      tau_at_lower_bound =
        abs(
          best_par["tau_days"] -
            tau_bounds[1L]
        ) <=
        tau_tolerance,
      tau_at_upper_bound =
        abs(
          best_par["tau_days"] -
            tau_bounds[2L]
        ) <=
        tau_tolerance,
      convergence = optimizer_convergence,
      optimizer_message = optimizer_message,
      matched_fit = matched_fit,
      prediction = final_prediction
    )
  )
}


#' Validate a permafrost SoilT model using leave-one-year-out validation
#'
#' Performs leave-one-year-out cross-validation for one model-site index and one
#' NEON soil-temperature depth. Each eligible calendar year is removed from the
#' likelihood, the requested process-model parameters are re-estimated from the
#' remaining years, and predictions are evaluated against observations in the
#' held-out year.
#'
#' When `fit_tau_only = TRUE`, every LOYO fold uses the same fixed values of
#' `a_C`, `n_warm`, and `n_cold`, and only `tau_days` is estimated from the
#' training years. This ensures that the held-out year does not influence tau.
#'
#' ERA5 forcing remains continuous across the complete requested time range for
#' every fold, so the process state is not reinitialized at calendar-year
#' boundaries.
#'
#' @md
#' @param lookup Site lookup table.
#' @param index Integer model-site index.
#' @param start_year First calendar year.
#' @param end_year Last calendar year.
#' @param multi_site Object returned by `get_soil_neon_data_multi_site()`.
#' @param era5_all Object returned by `get_soil_era5_data_multi_site()`.
#' @param vertical_position Character NEON vertical position.
#' @param fit_tau_only Logical. If `TRUE`, only tau is fitted in every fold.
#' @param fixed_a_C Fixed `a_C` in tau-only mode.
#' @param fixed_n_warm Fixed `n_warm` in tau-only mode.
#' @param fixed_n_cold Fixed `n_cold` in tau-only mode.
#' @param warmup_days Initial forcing spin-up excluded from validation.
#' @param min_observations Minimum matched 3-hour observations required for
#'   training and held-out evaluation.
#' @param min_days_per_year Minimum unique observation days required for an
#'   eligible held-out year.
#' @param min_days_per_season Minimum unique observation days required in each
#'   meteorological season of an eligible held-out year.
#' @param a_bounds Bounds for `a_C` in four-parameter mode.
#' @param n_warm_bounds Bounds for `n_warm` in four-parameter mode.
#' @param n_cold_bounds Bounds for `n_cold` in four-parameter mode.
#' @param tau_bounds Positive bounds for `tau_days`.
#' @param initial_parameters Optional starting parameter vector.
#' @param maxit Maximum L-BFGS-B iterations in four-parameter mode.
#'
#' @return A list containing year coverage, fold parameters, fold metrics,
#'   pooled out-of-sample metrics, and held-out predictions.
#'
#' @export
#' @author Yang Gu
validate_permafrost_loyo <- function(
    lookup,
    index,
    start_year,
    end_year,
    multi_site,
    era5_all,
    vertical_position = "501",
    fit_tau_only = TRUE,
    fixed_a_C = 0,
    fixed_n_warm = 1,
    fixed_n_cold = 0.3,
    warmup_days = 180L,
    min_observations = 100L,
    min_days_per_year = 120L,
    min_days_per_season = 10L,
    a_bounds = c(-10, 10),
    n_warm_bounds = c(0, 1.5),
    n_cold_bounds = c(0, 1.2),
    tau_bounds = c(0.125, 180),
    initial_parameters = NULL,
    maxit = 1000L
) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
  
  site_index <- as.integer(
    index
  )[1L]
  
  vertical_position <- as.character(
    vertical_position
  )[1L]
  
  start_year <- as.integer(
    start_year
  )[1L]
  
  end_year <- as.integer(
    end_year
  )[1L]

  warmup_days <- suppressWarnings(
    as.numeric(warmup_days)[1L]
  )

  min_observations <- suppressWarnings(
    as.integer(min_observations)[1L]
  )

  min_days_per_year <- suppressWarnings(
    as.integer(min_days_per_year)[1L]
  )

  min_days_per_season <- suppressWarnings(
    as.integer(min_days_per_season)[1L]
  )

  maxit <- suppressWarnings(
    as.integer(maxit)[1L]
  )

  if (!is.finite(site_index)) {
    stop("`index` must be one finite integer.", call. = FALSE)
  }

  if (
    !is.finite(start_year) ||
    !is.finite(end_year) ||
    end_year < start_year
  ) {
    stop("Invalid `start_year` / `end_year`.", call. = FALSE)
  }

  if (is.na(vertical_position) || !nzchar(vertical_position)) {
    stop("`vertical_position` must be one non-empty value.", call. = FALSE)
  }

  if (!is.finite(warmup_days) || warmup_days < 0) {
    stop("`warmup_days` must be one non-negative value.", call. = FALSE)
  }

  count_settings <- c(
    min_observations = min_observations,
    min_days_per_year = min_days_per_year,
    min_days_per_season = min_days_per_season,
    maxit = maxit
  )

  if (any(is.na(count_settings)) || any(count_settings < 1L)) {
    stop(
      "Observation/day thresholds and `maxit` must be positive integers.",
      call. = FALSE
    )
  }
  
  soilT_all <- getElement(
    multi_site,
    "soilT"
  )
  
  era5_mean_all <- getElement(
    era5_all,
    "era5_mean"
  )
  
  if (is.null(soilT_all)) {
    stop("`multi_site` does not contain `soilT`.", call. = FALSE)
  }
  
  if (is.null(era5_mean_all)) {
    stop("`era5_all` does not contain `era5_mean`.", call. = FALSE)
  }
  
  soilT_all <- data.table::as.data.table(
    data.table::copy(soilT_all)
  )
  
  era5_mean_all <- data.table::as.data.table(
    data.table::copy(era5_mean_all)
  )
  
  required_columns <- c(
    "time",
    "index",
    "verticalPosition",
    "SoilT"
  )
  
  if (!all(required_columns %in% names(soilT_all))) {
    stop(
      "`multi_site$soilT` is missing required columns.",
      call. = FALSE
    )
  }
  
  if (!all(c("time", "index") %in% names(era5_mean_all))) {
    stop(
      "`era5_all$era5_mean` must contain `time` and `index`.",
      call. = FALSE
    )
  }
  
  if (!"AirT_C" %in% names(era5_mean_all)) {
    if ("air_temperature" %in% names(era5_mean_all)) {
      era5_mean_all[
        ,
        AirT_C := as.numeric(air_temperature) - 273.15
      ]
    } else {
      stop(
        "`era5_all$era5_mean` must contain `AirT_C` or `air_temperature`.",
        call. = FALSE
      )
    }
  }
  
  soilT_all[
    ,
    time := as.POSIXct(time, tz = "UTC")
  ]
  
  era5_mean_all[
    ,
    time := as.POSIXct(time, tz = "UTC")
  ]
  
  start_time <- as.POSIXct(
    sprintf("%04d-01-01 00:00:00", start_year),
    tz = "UTC"
  )
  
  end_time <- as.POSIXct(
    sprintf("%04d-01-01 00:00:00", end_year + 1L),
    tz = "UTC"
  )
  
  obs_site <- soilT_all[
    index == site_index &
      as.character(verticalPosition) == vertical_position &
      time >= start_time &
      time < end_time &
      is.finite(SoilT)
  ]
  
  if (nrow(obs_site) == 0L) {
    stop(
      "No SoilT observations are available for LOYO validation.",
      call. = FALSE
    )
  }
  
  obs_site <- obs_site[
    ,
    .(
      Tsoil_obs_C = mean(SoilT, na.rm = TRUE)
    ),
    by = time
  ]
  
  met_times <- sort(
    unique(
      era5_mean_all[
        index == site_index &
          time >= start_time &
          time < end_time &
          !is.na(time) &
          is.finite(AirT_C),
        time
      ]
    )
  )
  
  if (length(met_times) < 2L) {
    stop(
      "Insufficient ERA5 forcing timestamps for LOYO validation.",
      call. = FALSE
    )
  }
  
  # Define coverage only from SoilT observations that can actually be matched
  # to usable ERA5 forcing, using the same warm-up origin as site fitting.
  warmup_end <-
    min(met_times) +
    as.numeric(warmup_days) *
    86400
  
  obs_site <- obs_site[
    time >= warmup_end &
      time %in% met_times
  ]
  
  if (nrow(obs_site) == 0L) {
    stop(
      "No matched SoilT observations remain after ERA5 matching and warm-up.",
      call. = FALSE
    )
  }
  
  obs_site[
    ,
    date :=
      as.Date(
        time,
        tz = "UTC"
      )
  ]
  
  obs_site[
    ,
    obs_year :=
      as.integer(
        format(
          time,
          "%Y",
          tz = "UTC"
        )
      )
  ]
  
  obs_site[
    ,
    month :=
      as.integer(
        format(
          time,
          "%m",
          tz = "UTC"
        )
      )
  ]
  
  obs_site[
    ,
    season :=
      ifelse(
        month %in%
          c(
            12L,
            1L,
            2L
          ),
        "DJF",
        ifelse(
          month %in%
            3L:5L,
          "MAM",
          ifelse(
            month %in%
              6L:8L,
            "JJA",
            "SON"
          )
        )
      )
  ]
  
  unique_days <- unique(
    obs_site[
      ,
      .(
        obs_year,
        date,
        season
      )
    ]
  )
  
  coverage <- unique_days[
    ,
    .(
      n_days =
        data.table::uniqueN(date),
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
    by = obs_year
  ]
  
  coverage[
    ,
    eligible :=
      n_days >=
      as.integer(min_days_per_year) &
      n_DJF >=
      as.integer(min_days_per_season) &
      n_MAM >=
      as.integer(min_days_per_season) &
      n_JJA >=
      as.integer(min_days_per_season) &
      n_SON >=
      as.integer(min_days_per_season)
  ]
  
  data.table::setorder(
    coverage,
    obs_year
  )
  
  eligible_years <- coverage[
    eligible == TRUE,
    obs_year
  ]
  
  if (length(eligible_years) < 2L) {
    stop(
      "LOYO requires at least two eligible years. Eligible years: ",
      paste(
        eligible_years,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  all_obs_years <- sort(
    unique(
      obs_site$obs_year
    )
  )
  
  fold_parameter_rows <- list()
  fold_metric_rows <- list()
  validation_rows <- list()
  
  for (held_out_year in eligible_years) {
    
    train_years <- all_obs_years[
      all_obs_years !=
        held_out_year
    ]
    
    fold_fit <- fit_permafrost_site_mle(
      lookup = lookup,
      index = site_index,
      start_year = start_year,
      end_year = end_year,
      multi_site = multi_site,
      era5_all = era5_all,
      vertical_position = vertical_position,
      fit_years = train_years,
      fit_tau_only = fit_tau_only,
      fixed_a_C = fixed_a_C,
      fixed_n_warm = fixed_n_warm,
      fixed_n_cold = fixed_n_cold,
      warmup_days = warmup_days,
      min_observations = min_observations,
      a_bounds = a_bounds,
      n_warm_bounds = n_warm_bounds,
      n_cold_bounds = n_cold_bounds,
      tau_bounds = tau_bounds,
      initial_parameters = initial_parameters,
      maxit = maxit
    )
    
    fold_parameters <- getElement(
      fold_fit,
      "parameters"
    )
    
    prediction <- getElement(
      fold_fit,
      "prediction"
    )
    
    test_obs <- obs_site[
      obs_year == held_out_year,
      .(
        time,
        Tsoil_obs_C
      )
    ]
    
    prediction_row <- match(
      test_obs$time,
      prediction$time
    )
    
    test_obs[
      ,
      prediction_row :=
        prediction_row
    ]
    
    test_obs <- test_obs[
      !is.na(prediction_row)
    ]
    
    if (nrow(test_obs) == 0L) {
      next
    }
    
    prediction_values <- getElement(
      prediction,
      "Tsoil_pred_C"
    )
    
    test_prediction <- prediction_values[
      test_obs$prediction_row
    ]
    
    test_obs[
      ,
      Tsoil_pred_C :=
        test_prediction
    ]
    
    fold_metrics <- permafrost_tsoil_data_model(
      observed_C = test_obs$Tsoil_obs_C,
      predicted_C = test_obs$Tsoil_pred_C,
      min_observations = min_observations
    )
    
    if (!is.finite(getElement(fold_metrics, "nll"))) {
      next
    }
    
    fold_parameter_rows[
      length(fold_parameter_rows) +
        1L
    ] <- list(
      data.table::data.table(
        held_out_year = held_out_year,
        fit_mode = getElement(fold_fit, "fit_mode"),
        a_C = as.numeric(fold_parameters["a_C"]),
        n_warm = as.numeric(fold_parameters["n_warm"]),
        n_cold = as.numeric(fold_parameters["n_cold"]),
        tau_days = as.numeric(fold_parameters["tau_days"]),
        tau_at_lower_bound = getElement(fold_fit, "tau_at_lower_bound"),
        tau_at_upper_bound = getElement(fold_fit, "tau_at_upper_bound"),
        n_train = getElement(fold_fit, "n_fit"),
        training_R2 = getElement(fold_fit, "training_r2")
      )
    )
    
    fold_metric_rows[
      length(fold_metric_rows) +
        1L
    ] <- list(
      data.table::data.table(
        held_out_year = held_out_year,
        n = getElement(fold_metrics, "n"),
        r2 = getElement(fold_metrics, "r2"),
        rmse_C = getElement(fold_metrics, "rmse_C"),
        mae_C = getElement(fold_metrics, "mae_C"),
        bias_C = getElement(fold_metrics, "bias_C"),
        correlation = getElement(fold_metrics, "correlation")
      )
    )
    
    test_obs[
      ,
      `:=`(
        held_out_year = held_out_year,
        fold_a_C = as.numeric(fold_parameters["a_C"]),
        fold_n_warm = as.numeric(fold_parameters["n_warm"]),
        fold_n_cold = as.numeric(fold_parameters["n_cold"]),
        fold_tau_days = as.numeric(fold_parameters["tau_days"])
      )
    ]
    
    test_obs[
      ,
      prediction_row :=
        NULL
    ]
    
    validation_rows[
      length(validation_rows) +
        1L
    ] <- list(
      test_obs
    )
  }
  
  if (length(validation_rows) < 2L) {
    stop(
      "LOYO requires at least two valid held-out folds after ERA5 matching ",
      "and metric checks; only ",
      length(validation_rows),
      " valid fold(s) were produced.",
      call. = FALSE
    )
  }
  
  fold_parameters <- data.table::rbindlist(
    fold_parameter_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  fold_metrics <- data.table::rbindlist(
    fold_metric_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  validation <- data.table::rbindlist(
    validation_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  data.table::setorder(
    validation,
    time
  )
  
  pooled_metrics_raw <- permafrost_tsoil_data_model(
    observed_C = validation$Tsoil_obs_C,
    predicted_C = validation$Tsoil_pred_C,
    min_observations = min_observations
  )
  
  pooled_metrics <- data.table::data.table(
    n = getElement(pooled_metrics_raw, "n"),
    n_folds = data.table::uniqueN(validation$held_out_year),
    r2 = getElement(pooled_metrics_raw, "r2"),
    rmse_C = getElement(pooled_metrics_raw, "rmse_C"),
    mae_C = getElement(pooled_metrics_raw, "mae_C"),
    bias_C = getElement(pooled_metrics_raw, "bias_C"),
    correlation = getElement(pooled_metrics_raw, "correlation")
  )
  
  return(
    list(
      year_coverage = coverage,
      eligible_years = eligible_years,
      fold_parameters = fold_parameters,
      fold_metrics = fold_metrics,
      pooled_metrics = pooled_metrics,
      validation = validation
    )
  )
}


#' Fit and validate the permafrost SoilT model for all requested sites
#'
#' Runs maximum-likelihood fitting and leave-one-year-out validation for every
#' requested model-site index and NEON soil-temperature vertical position found
#' in pre-extracted calibration data.
#'
#' The `lookup` object controls which site indices are processed. The function
#' does not independently classify or detect permafrost; therefore, when the
#' supplied lookup contains only permafrost sites, every processed job is a
#' permafrost-site calibration.
#'
#' With `fit_tau_only = TRUE`, the same fixed `a_C`, `n_warm`, and `n_cold`
#' values are used for every site and every LOYO fold, and only site-specific
#' `tau_days` values are estimated. With `fit_tau_only = FALSE`, all four
#' process parameters are fitted independently for every site and LOYO fold.
#'
#' @md
#' @param lookup Site lookup table containing the requested `index` values.
#' @param start_year First calendar year.
#' @param end_year Last calendar year.
#' @param multi_site Object returned by `get_soil_neon_data_multi_site()`.
#' @param era5_all Object returned by `get_soil_era5_data_multi_site()`.
#' @param vertical_positions Character vector of requested NEON vertical
#'   positions. `NULL` uses every vertical position available for the requested
#'   lookup indices in `multi_site$soilT`.
#' @param workers Number of parallel workers. Set to one for sequential fitting.
#' @param fit_tau_only Logical. If `TRUE`, estimate tau only.
#' @param fixed_a_C Fixed baseline soil-air offset in tau-only mode.
#' @param fixed_n_warm Fixed warm coupling coefficient in tau-only mode.
#' @param fixed_n_cold Fixed cold coupling coefficient in tau-only mode.
#' @param warmup_days Initial forcing spin-up excluded from likelihood and LOYO.
#' @param min_observations Minimum matched 3-hour observations required.
#' @param min_days_per_year Minimum unique observation days for an eligible
#'   LOYO held-out year.
#' @param min_days_per_season Minimum unique observation days in every season
#'   for an eligible LOYO held-out year.
#' @param a_bounds Bounds for `a_C` in four-parameter mode.
#' @param n_warm_bounds Bounds for `n_warm` in four-parameter mode.
#' @param n_cold_bounds Bounds for `n_cold` in four-parameter mode.
#' @param tau_bounds Positive lower and upper bounds for tau in days.
#' @param initial_parameters Optional starting parameter vector.
#' @param maxit Maximum L-BFGS-B iterations in four-parameter mode.
#' @param output_dir Optional directory for combined output CSV files. Set to
#'   `NULL` to suppress file writing.
#'
#' @return A list containing combined final-fit summary, LOYO year coverage,
#'   fold parameters, fold metrics, held-out predictions, and processing status.
#'
#' @export
#' @author Yang Gu
fit_permafrost_all_sites <- function(
    lookup,
    start_year,
    end_year,
    multi_site,
    era5_all,
    vertical_positions = NULL,
    workers = 1L,
    fit_tau_only = TRUE,
    fixed_a_C = 0,
    fixed_n_warm = 1,
    fixed_n_cold = 0.3,
    warmup_days = 180L,
    min_observations = 100L,
    min_days_per_year = 120L,
    min_days_per_season = 10L,
    a_bounds = c(-10, 10),
    n_warm_bounds = c(0, 1.5),
    n_cold_bounds = c(0, 1.2),
    tau_bounds = c(0.125, 180),
    initial_parameters = NULL,
    maxit = 1000L,
    output_dir = NULL
) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }
  
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
  
  if (!"index" %in% names(lookup_dt)) {
    stop("`lookup` must contain `index`.", call. = FALSE)
  }
  
  target_indices <- sort(
    unique(
      as.integer(
        lookup_dt[
          !is.na(index),
          index
        ]
      )
    )
  )
  
  target_indices <- target_indices[
    is.finite(target_indices)
  ]
  
  if (length(target_indices) == 0L) {
    stop("No valid indices were found in `lookup`.", call. = FALSE)
  }
  
  start_year <- suppressWarnings(as.integer(start_year)[1L])
  end_year <- suppressWarnings(as.integer(end_year)[1L])
  
  if (
    !is.finite(start_year) ||
    !is.finite(end_year) ||
    end_year < start_year
  ) {
    stop("Invalid `start_year` / `end_year`.", call. = FALSE)
  }
  
  if (
    length(warmup_days) != 1L ||
    !is.finite(warmup_days) ||
    warmup_days < 0
  ) {
    stop("`warmup_days` must be one non-negative finite value.", call. = FALSE)
  }
  
  count_settings <- c(
    min_observations = min_observations,
    min_days_per_year = min_days_per_year,
    min_days_per_season = min_days_per_season
  )
  
  if (any(!is.finite(count_settings)) || any(count_settings < 1)) {
    stop(
      "Observation/day threshold arguments must be positive finite values.",
      call. = FALSE
    )
  }
  
  soilT_all <- getElement(
    multi_site,
    "soilT"
  )
  
  era5_mean_all <- getElement(
    era5_all,
    "era5_mean"
  )
  
  if (is.null(soilT_all)) {
    stop("`multi_site` does not contain `soilT`.", call. = FALSE)
  }
  
  if (is.null(era5_mean_all)) {
    stop("`era5_all` does not contain `era5_mean`.", call. = FALSE)
  }
  
  soilT_all <- data.table::as.data.table(
    data.table::copy(
      soilT_all
    )
  )
  
  era5_mean_all <- data.table::as.data.table(
    data.table::copy(
      era5_mean_all
    )
  )
  
  required_soil_columns <- c(
    "index",
    "verticalPosition",
    "time",
    "SoilT"
  )
  
  if (
    !all(
      required_soil_columns %in%
      names(soilT_all)
    )
  ) {
    stop(
      "`multi_site$soilT` is missing required columns.",
      call. = FALSE
    )
  }
  
  if (
    !all(
      c(
        "index",
        "time"
      ) %in%
      names(era5_mean_all)
    )
  ) {
    stop(
      "`era5_all$era5_mean` must contain `index` and `time`.",
      call. = FALSE
    )
  }
  
  soilT_all[
    ,
    time := as.POSIXct(time, tz = "UTC")
  ]
  
  start_time <- as.POSIXct(
    sprintf("%04d-01-01 00:00:00", start_year),
    tz = "UTC"
  )
  
  end_time <- as.POSIXct(
    sprintf("%04d-01-01 00:00:00", end_year + 1L),
    tz = "UTC"
  )
  
  soilT_target <- soilT_all[
    index %in%
      target_indices &
      time >= start_time &
      time < end_time &
      is.finite(SoilT)
  ]
  
  if (!is.null(vertical_positions)) {
    
    vertical_positions <- as.character(
      vertical_positions
    )
    
    soilT_target <- soilT_target[
      as.character(verticalPosition) %in%
        vertical_positions
    ]
  }
  
  jobs <- unique(
    soilT_target[
      ,
      .(
        index,
        verticalPosition =
          as.character(verticalPosition)
      )
    ]
  )
  
  data.table::setorder(
    jobs,
    index,
    verticalPosition
  )
  
  if (nrow(jobs) == 0L) {
    stop(
      paste0(
        "No site-depth combinations in `multi_site$soilT` match ",
        "the requested lookup and vertical positions."
      ),
      call. = FALSE
    )
  }
  
  fit_mode <- if (isTRUE(fit_tau_only)) {
    "tau_only"
  } else {
    "four_parameter"
  }
  
  message(
    "SoilT MLE: ",
    nrow(jobs),
    " site-depth jobs across ",
    length(target_indices),
    " lookup indices | mode = ",
    fit_mode,
    "."
  )
  
  if (isTRUE(fit_tau_only)) {
    message(
      "Fixed parameters: a_C = ",
      fixed_a_C,
      ", n_warm = ",
      fixed_n_warm,
      ", n_cold = ",
      fixed_n_cold,
      "."
    )
  }
  
  run_one_job <- function(
    job_number
  ) {
    
    site_index <- jobs$index[job_number]
    vertical_position <- jobs$verticalPosition[job_number]
    
    site_lookup <- lookup_dt[
      index == site_index
    ]
    
    local_multi_site <- list(
      soilT = soilT_target[
        index == site_index
      ]
    )
    
    local_era5_all <- list(
      era5_mean = era5_mean_all[
        index == site_index
      ]
    )
    
    base_args <- list(
      lookup = site_lookup,
      index = site_index,
      start_year = start_year,
      end_year = end_year,
      multi_site = local_multi_site,
      era5_all = local_era5_all,
      vertical_position = vertical_position,
      fit_tau_only = fit_tau_only,
      fixed_a_C = fixed_a_C,
      fixed_n_warm = fixed_n_warm,
      fixed_n_cold = fixed_n_cold,
      warmup_days = warmup_days,
      min_observations = min_observations,
      a_bounds = a_bounds,
      n_warm_bounds = n_warm_bounds,
      n_cold_bounds = n_cold_bounds,
      tau_bounds = tau_bounds,
      initial_parameters = initial_parameters,
      maxit = maxit
    )
    
    # Estimate the final all-year parameters independently of LOYO. A site can
    # therefore retain a valid final MLE even when temporal validation fails.
    final_fit <- tryCatch(
      do.call(
        fit_permafrost_site_mle,
        c(
          base_args,
          list(fit_years = NULL)
        )
      ),
      error = function(e) e
    )
    
    if (inherits(final_fit, "error")) {
      final_message <- conditionMessage(final_fit)
      message(
        "FINAL FIT FAILED | index=",
        site_index,
        " | VER=",
        vertical_position,
        " | ",
        final_message
      )
      
      return(
        list(
          success = FALSE,
          final_fit_success = FALSE,
          loyo_success = FALSE,
          summary = NULL,
          year_coverage = NULL,
          fold_parameters = NULL,
          fold_metrics = NULL,
          validation = NULL,
          message = final_message,
          loyo_message = NA_character_
        )
      )
    }
    
    loyo <- tryCatch(
      do.call(
        validate_permafrost_loyo,
        c(
          base_args,
          list(
            min_days_per_year = min_days_per_year,
            min_days_per_season = min_days_per_season
          )
        )
      ),
      error = function(e) e
    )
    
    loyo_success <- !inherits(loyo, "error")
    loyo_message <- if (loyo_success) {
      "OK"
    } else {
      conditionMessage(loyo)
    }
    
    if (!loyo_success) {
      message(
        "LOYO FAILED; FINAL FIT RETAINED | index=",
        site_index,
        " | VER=",
        vertical_position,
        " | ",
        loyo_message
      )
    }
    
    final_parameters <- getElement(
      final_fit,
      "parameters"
    )
    
    pooled_metrics <- if (loyo_success) {
      getElement(loyo, "pooled_metrics")
    } else {
      data.table::data.table(
        n = 0L,
        n_folds = 0L,
        r2 = NA_real_,
        rmse_C = NA_real_,
        mae_C = NA_real_,
        bias_C = NA_real_,
        correlation = NA_real_
      )
    }

    fold_tau_values <- if (loyo_success) {
      as.numeric(
        getElement(
          getElement(loyo, "fold_parameters"),
          "tau_days"
        )
      )
    } else {
      numeric()
    }

    fold_tau_values <- fold_tau_values[
      is.finite(fold_tau_values)
    ]

    fold_tau_mean <- if (length(fold_tau_values) > 0L) {
      mean(fold_tau_values)
    } else {
      NA_real_
    }

    fold_tau_sd <- if (length(fold_tau_values) > 1L) {
      stats::sd(fold_tau_values)
    } else {
      NA_real_
    }

    fold_tau_cv <- if (
      is.finite(fold_tau_mean) &&
      fold_tau_mean != 0 &&
      is.finite(fold_tau_sd)
    ) {
      fold_tau_sd / abs(fold_tau_mean)
    } else {
      NA_real_
    }
    
    summary <- data.table::data.table(
      index = site_index,
      AmeriFlux_ID = getElement(final_fit, "AmeriFlux_ID"),
      NEON_code = getElement(final_fit, "NEON_code"),
      verticalPosition = vertical_position,
      mean_zOffset = getElement(final_fit, "mean_zOffset"),
      fit_mode = getElement(final_fit, "fit_mode"),
      fixed_a_C = getElement(final_fit, "fixed_a_C"),
      fixed_n_warm = getElement(final_fit, "fixed_n_warm"),
      fixed_n_cold = getElement(final_fit, "fixed_n_cold"),
      a_C = as.numeric(final_parameters["a_C"]),
      n_warm = as.numeric(final_parameters["n_warm"]),
      n_cold = as.numeric(final_parameters["n_cold"]),
      tau_days = as.numeric(final_parameters["tau_days"]),
      tau_at_lower_bound = getElement(final_fit, "tau_at_lower_bound"),
      tau_at_upper_bound = getElement(final_fit, "tau_at_upper_bound"),
      residual_sigma_C = getElement(final_fit, "sigma_C"),
      final_nll = getElement(final_fit, "nll"),
      n_final_fit = getElement(final_fit, "n_fit"),
      training_R2 = getElement(final_fit, "training_r2"),
      training_RMSE_C = getElement(final_fit, "training_rmse_C"),
      optimizer_convergence = getElement(final_fit, "convergence"),
      optimizer_message = getElement(final_fit, "optimizer_message"),
      LOYO_success = loyo_success,
      LOYO_message = loyo_message,
      LOYO_n = getElement(pooled_metrics, "n")[1L],
      LOYO_n_folds = getElement(pooled_metrics, "n_folds")[1L],
      LOYO_tau_mean_days = fold_tau_mean,
      LOYO_tau_SD_days = fold_tau_sd,
      LOYO_tau_CV = fold_tau_cv,
      LOYO_tau_min_days = if (length(fold_tau_values) > 0L) {
        min(fold_tau_values)
      } else {
        NA_real_
      },
      LOYO_tau_max_days = if (length(fold_tau_values) > 0L) {
        max(fold_tau_values)
      } else {
        NA_real_
      },
      LOYO_tau_range_days = if (length(fold_tau_values) > 0L) {
        diff(range(fold_tau_values))
      } else {
        NA_real_
      },
      LOYO_R2 = getElement(pooled_metrics, "r2")[1L],
      LOYO_RMSE_C = getElement(pooled_metrics, "rmse_C")[1L],
      LOYO_MAE_C = getElement(pooled_metrics, "mae_C")[1L],
      LOYO_bias_C = getElement(pooled_metrics, "bias_C")[1L],
      LOYO_correlation = getElement(pooled_metrics, "correlation")[1L]
    )
    
    fold_parameters <- if (loyo_success) {
      data.table::copy(getElement(loyo, "fold_parameters"))
    } else {
      data.table::data.table()
    }
    
    fold_metrics <- if (loyo_success) {
      data.table::copy(getElement(loyo, "fold_metrics"))
    } else {
      data.table::data.table()
    }
    
    validation <- if (loyo_success) {
      data.table::copy(getElement(loyo, "validation"))
    } else {
      data.table::data.table()
    }
    
    year_coverage <- if (loyo_success) {
      data.table::copy(getElement(loyo, "year_coverage"))
    } else {
      data.table::data.table()
    }
    
    if (nrow(year_coverage) > 0L) {
      year_coverage[
        ,
        `:=`(
          index = site_index,
          verticalPosition = vertical_position
        )
      ]
    }
    
    if (nrow(fold_parameters) > 0L) {
      fold_parameters[
        ,
        `:=`(
          index = site_index,
          verticalPosition = vertical_position
        )
      ]
    }
    
    if (nrow(fold_metrics) > 0L) {
      fold_metrics[
        ,
        `:=`(
          index = site_index,
          verticalPosition = vertical_position
        )
      ]
    }
    
    if (nrow(validation) > 0L) {
      validation[
        ,
        `:=`(
          index = site_index,
          verticalPosition = vertical_position
        )
      ]
    }
    
    list(
      success = TRUE,
      final_fit_success = TRUE,
      loyo_success = loyo_success,
      summary = summary,
      year_coverage = year_coverage,
      fold_parameters = fold_parameters,
      fold_metrics = fold_metrics,
      validation = validation,
      message = if (loyo_success) {
        "OK"
      } else {
        paste0("Final fit succeeded; LOYO failed: ", loyo_message)
      },
      loyo_message = loyo_message
    )
  }

  workers <- as.integer(
    workers
  )[1L]
  
  if (
    !is.finite(workers) ||
    workers < 1L
  ) {
    stop("`workers` must be >= 1.", call. = FALSE)
  }
  
  job_numbers <- seq_len(
    nrow(jobs)
  )
  
  if (workers == 1L) {
    
    job_results <- lapply(
      job_numbers,
      run_one_job
    )
    
  } else {
    
    if (!requireNamespace("future", quietly = TRUE)) {
      stop(
        "Package `future` is required for parallel execution.",
        call. = FALSE
      )
    }
    
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop(
        "Package `furrr` is required for parallel execution.",
        call. = FALSE
      )
    }
    
    previous_future_plan <- future::plan()
    
    on.exit(
      future::plan(previous_future_plan),
      add = TRUE
    )
    
    future::plan(
      future::multisession,
      workers = workers
    )
    
    job_results <- furrr::future_map(
      job_numbers,
      run_one_job,
      .options = furrr::furrr_options(
        seed = TRUE
      )
    )
  }
  
  summary_rows <- list()
  year_coverage_rows <- list()
  fold_parameter_rows <- list()
  fold_metric_rows <- list()
  validation_rows <- list()
  status_rows <- list()
  
  for (job_i in seq_along(job_results)) {
    
    result_i <- do.call(
      c,
      job_results[job_i]
    )
    
    success_i <- isTRUE(
      getElement(
        result_i,
        "success"
      )
    )
    
    final_fit_success_i <- isTRUE(
      getElement(result_i, "final_fit_success")
    )
    
    loyo_success_i <- isTRUE(
      getElement(result_i, "loyo_success")
    )
    
    status_rows[
      length(status_rows) +
        1L
    ] <- list(
      data.table::data.table(
        index = jobs$index[job_i],
        verticalPosition = jobs$verticalPosition[job_i],
        success = success_i,
        final_fit_success = final_fit_success_i,
        loyo_success = loyo_success_i,
        message = getElement(result_i, "message"),
        loyo_message = getElement(result_i, "loyo_message")
      )
    )
    
    if (!success_i) {
      next
    }
    
    summary_rows[
      length(summary_rows) +
        1L
    ] <- list(
      getElement(result_i, "summary")
    )
    
    year_coverage_i <- getElement(result_i, "year_coverage")
    if (!is.null(year_coverage_i) && nrow(year_coverage_i) > 0L) {
      year_coverage_rows[
        length(year_coverage_rows) + 1L
      ] <- list(year_coverage_i)
    }
    
    fold_parameter_rows[
      length(fold_parameter_rows) +
        1L
    ] <- list(
      getElement(result_i, "fold_parameters")
    )
    
    fold_metric_rows[
      length(fold_metric_rows) +
        1L
    ] <- list(
      getElement(result_i, "fold_metrics")
    )
    
    validation_rows[
      length(validation_rows) +
        1L
    ] <- list(
      getElement(result_i, "validation")
    )
  }
  
  summary_table <- if (length(summary_rows) > 0L) {
    data.table::rbindlist(
      summary_rows,
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table::data.table()
  }
  
  year_coverage <- if (length(year_coverage_rows) > 0L) {
    data.table::rbindlist(
      year_coverage_rows,
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table::data.table()
  }
  
  fold_parameters <- if (length(fold_parameter_rows) > 0L) {
    data.table::rbindlist(
      fold_parameter_rows,
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table::data.table()
  }
  
  fold_metrics <- if (length(fold_metric_rows) > 0L) {
    data.table::rbindlist(
      fold_metric_rows,
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table::data.table()
  }
  
  validation <- if (length(validation_rows) > 0L) {
    data.table::rbindlist(
      validation_rows,
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table::data.table()
  }
  
  status <- data.table::rbindlist(
    status_rows,
    use.names = TRUE,
    fill = TRUE
  )
  
  if (nrow(summary_table) > 0L) {
    data.table::setorder(
      summary_table,
      index,
      verticalPosition
    )
  }

  if (nrow(year_coverage) > 0L) {
    data.table::setorder(
      year_coverage,
      index,
      verticalPosition,
      obs_year
    )
  }

  if (nrow(fold_parameters) > 0L) {
    data.table::setorder(
      fold_parameters,
      index,
      verticalPosition,
      held_out_year
    )
  }

  if (nrow(fold_metrics) > 0L) {
    data.table::setorder(
      fold_metrics,
      index,
      verticalPosition,
      held_out_year
    )
  }

  if (nrow(validation) > 0L) {
    data.table::setorder(
      validation,
      index,
      verticalPosition,
      held_out_year,
      time
    )
  }
  
  if (nrow(status) > 0L) {
    data.table::setorder(
      status,
      index,
      verticalPosition
    )
  }
  
  if (!is.null(output_dir)) {
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    mode_tag <- if (isTRUE(fit_tau_only)) {
      "tau_only"
    } else {
      "four_parameter"
    }
    
    directory_label <- basename(
      normalizePath(output_dir, mustWork = FALSE)
    )
    
    output_tag <- if (
      directory_label %in% c("non_permafrost", "permafrost")
    ) {
      directory_label
    } else {
      "tsoil"
    }
    
    if (ncol(summary_table) > 0L) {
      data.table::fwrite(
        summary_table,
        file.path(
          output_dir,
          paste0(
            output_tag,
            "_tsoil_",
            mode_tag,
            "_MLE_summary.csv"
          )
        )
      )
    }
    
    if (ncol(year_coverage) > 0L) {
      data.table::fwrite(
        year_coverage,
        file.path(
          output_dir,
          paste0(
            output_tag,
            "_tsoil_",
            mode_tag,
            "_LOYO_year_coverage.csv"
          )
        )
      )
    }
    
    if (ncol(fold_parameters) > 0L) {
      data.table::fwrite(
        fold_parameters,
        file.path(
          output_dir,
          paste0(
            output_tag,
            "_tsoil_",
            mode_tag,
            "_LOYO_fold_parameters.csv"
          )
        )
      )
    }
    
    if (ncol(fold_metrics) > 0L) {
      data.table::fwrite(
        fold_metrics,
        file.path(
          output_dir,
          paste0(
            output_tag,
            "_tsoil_",
            mode_tag,
            "_LOYO_fold_metrics.csv"
          )
        )
      )
    }
    
    if (ncol(validation) > 0L) {
      data.table::fwrite(
        validation,
        file.path(
          output_dir,
          paste0(
            output_tag,
            "_tsoil_",
            mode_tag,
            "_LOYO_predictions.csv"
          )
        )
      )
    }
    
    data.table::fwrite(
      status,
      file.path(
        output_dir,
        paste0(
          output_tag,
          "_tsoil_",
          mode_tag,
          "_status.csv"
        )
      )
    )
  }
  
  message(
    "SoilT MLE finished: ",
    sum(status$success),
    " successful / ",
    nrow(status),
    " requested site-depth jobs | mode = ",
    fit_mode,
    "."
  )
  
  return(
    list(
      summary = summary_table,
      year_coverage = year_coverage,
      fold_parameters = fold_parameters,
      fold_metrics = fold_metrics,
      validation = validation,
      status = status
    )
  )
}


# =============================================================================
# Example: tau-only permafrost SoilT calibration
#
# permafrost_tau_results <- fit_permafrost_all_sites(
#   lookup = permafrost,
#   start_year = 2017L,
#   end_year = 2024L,
#   multi_site = multi_site,
#   era5_all = era5_all,
#   vertical_positions = "501",
#   workers = 5L,
#   fit_tau_only = TRUE,
#   fixed_a_C = 0,
#   fixed_n_warm = 1,
#   fixed_n_cold = 0.3,
#   tau_bounds = c(0.125, 180),
#   warmup_days = 180L,
#   min_observations = 100L,
#   min_days_per_year = 120L,
#   min_days_per_season = 10L,
#   output_dir = "/projectnb/dietzelab/guYANG/soilparam/permafrost_tau_only"
# )
# =============================================================================


# =============================================================================
# Unified NEON soil-temperature MLE workflow
# =============================================================================


#' Fit site-specific soil-temperature tau for all NEON sites
#'
#' Splits a NEON site lookup into permafrost and non-permafrost sites using
#' latitude, fits each site/depth independently, and performs leave-one-year-out
#' (LOYO) validation. Sites above `permafrost_latitude_threshold` use the
#' asymmetric permafrost process model. All other sites use the standard
#' non-permafrost tau model, represented by the same process equation with
#' `a_C = 0`, `n_warm = 1`, and `n_cold = 1`, so effective air temperature is
#' exactly observed air temperature.
#'
#' The function delegates fitting and LOYO validation to
#' [fit_permafrost_all_sites()]. Despite that historical function name, fixing
#' the three air-coupling parameters to `0, 1, 1` gives the non-permafrost
#' causal tau model exactly. Permafrost parameters may either be fixed while
#' fitting tau only or estimated jointly with tau.
#'
#' A latitude exactly equal to the threshold is classified as non-permafrost.
#' Missing or non-finite latitudes are rejected rather than silently assigned
#' to a process model.
#'
#' @param lookup A data.frame or data.table with one row per site and columns
#'   `index` and `site_lat`. Other site metadata are retained by the underlying
#'   fitting workflow.
#' @param start_year,end_year First and last calendar years used for fitting.
#' @param multi_site Object returned by `get_soil_neon_data_multi_site()`; it
#'   must contain `soilT`.
#' @param era5_all Object returned by `get_soil_era5_data_multi_site()`; it must
#'   contain `era5_mean`.
#' @param vertical_positions Character vector of NEON soil-temperature depth
#'   codes. `NULL` uses every available depth.
#' @param workers Number of parallel workers used separately for each process
#'   model group.
#' @param permafrost_latitude_threshold Sites with `site_lat` strictly greater
#'   than this value use the permafrost process model. Default is 63 degrees.
#' @param permafrost_fit_tau_only If `TRUE`, only permafrost tau is fitted.
#' @param permafrost_fixed_a_C,permafrost_fixed_n_warm,permafrost_fixed_n_cold
#'   Fixed permafrost parameters used in tau-only mode.
#' @param warmup_days Initial forcing spin-up excluded from fitting and LOYO.
#' @param min_observations Minimum matched observations required.
#' @param min_days_per_year Minimum unique observation days required for an
#'   eligible held-out year.
#' @param min_days_per_season Minimum unique observation days required in every
#'   meteorological season of an eligible held-out year.
#' @param a_bounds,n_warm_bounds,n_cold_bounds Bounds used when jointly fitting
#'   all permafrost process parameters.
#' @param tau_bounds Positive lower and upper bounds for tau in days, used for
#'   both process-model groups.
#' @param initial_parameters Optional initial parameters for a joint
#'   permafrost fit.
#' @param maxit Maximum optimizer iterations.
#' @param output_dir Optional output directory. When supplied, group-specific
#'   detailed files are written under `non_permafrost/` and `permafrost/`, and
#'   the combined best-tau table is written at the top level.
#'
#' @return A list with exactly three top-level elements:
#'   * `best_tau`: successful final-fit rows for every site/depth, including
#'     `process_model` and `is_permafrost` classifications.
#'   * `non_permafrost_loyo`: LOYO fold parameters, metrics, predictions, and
#'     status for non-permafrost sites.
#'   * `permafrost_loyo`: the corresponding LOYO results for permafrost sites.
#'
#' @export
#' @author Yang Gu
neon_tsoil_mle <- function(
    lookup,
    start_year,
    end_year,
    multi_site,
    era5_all,
    vertical_positions = NULL,
    workers = 1L,
    permafrost_latitude_threshold = 63,
    permafrost_fit_tau_only = TRUE,
    permafrost_fixed_a_C = 0,
    permafrost_fixed_n_warm = 1,
    permafrost_fixed_n_cold = 0.3,
    warmup_days = 180L,
    min_observations = 100L,
    min_days_per_year = 120L,
    min_days_per_season = 10L,
    a_bounds = c(-10, 10),
    n_warm_bounds = c(0, 1.5),
    n_cold_bounds = c(0, 1.2),
    tau_bounds = c(0.125, 180),
    initial_parameters = NULL,
    maxit = 1000L,
    output_dir = NULL
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required.", call. = FALSE)
  }

  if (!exists(
    "fit_permafrost_all_sites",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      paste0(
        "`fit_permafrost_all_sites()` is required. Source ",
        "this standalone workflow file completely before calling ",
        "`neon_tsoil_mle()`."
      ),
      call. = FALSE
    )
  }

  lookup_dt <- data.table::as.data.table(
    data.table::copy(lookup)
  )

  required_lookup_columns <- c("index", "site_lat")
  missing_lookup_columns <- setdiff(
    required_lookup_columns,
    names(lookup_dt)
  )

  if (length(missing_lookup_columns) > 0L) {
    stop(
      "`lookup` is missing required columns: ",
      paste(missing_lookup_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  lookup_dt[
    ,
    `:=`(
      index = suppressWarnings(as.integer(index)),
      site_lat = suppressWarnings(as.numeric(site_lat))
    )
  ]

  if (any(!is.finite(lookup_dt$index))) {
    stop("`lookup$index` contains missing or invalid values.", call. = FALSE)
  }

  if (anyDuplicated(lookup_dt$index) > 0L) {
    stop("`lookup` must contain exactly one row per `index`.", call. = FALSE)
  }

  if (any(!is.finite(lookup_dt$site_lat))) {
    invalid_indices <- lookup_dt[
      !is.finite(site_lat),
      index
    ]
    stop(
      "`site_lat` is missing or invalid for index: ",
      paste(invalid_indices, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  latitude_threshold <- suppressWarnings(
    as.numeric(permafrost_latitude_threshold)[1L]
  )
  if (!is.finite(latitude_threshold)) {
    stop(
      "`permafrost_latitude_threshold` must be finite.",
      call. = FALSE
    )
  }

  lookup_dt[
    ,
    `:=`(
      is_permafrost = site_lat > latitude_threshold,
      process_model = ifelse(
        site_lat > latitude_threshold,
        "permafrost",
        "non_permafrost"
      )
    )
  ]

  empty_loyo_result <- function() {
    list(
      year_coverage = data.table::data.table(),
      fold_parameters = data.table::data.table(),
      fold_metrics = data.table::data.table(),
      validation = data.table::data.table(),
      status = data.table::data.table()
    )
  }

  run_group <- function(group_lookup, process_model, group_output_dir) {
    if (nrow(group_lookup) == 0L) {
      return(NULL)
    }

    is_permafrost_group <- identical(process_model, "permafrost")

    fit_permafrost_all_sites(
      lookup = group_lookup,
      start_year = start_year,
      end_year = end_year,
      multi_site = multi_site,
      era5_all = era5_all,
      vertical_positions = vertical_positions,
      workers = workers,
      fit_tau_only = if (is_permafrost_group) {
        permafrost_fit_tau_only
      } else {
        TRUE
      },
      fixed_a_C = if (is_permafrost_group) {
        permafrost_fixed_a_C
      } else {
        0
      },
      fixed_n_warm = if (is_permafrost_group) {
        permafrost_fixed_n_warm
      } else {
        1
      },
      fixed_n_cold = if (is_permafrost_group) {
        permafrost_fixed_n_cold
      } else {
        1
      },
      warmup_days = warmup_days,
      min_observations = min_observations,
      min_days_per_year = min_days_per_year,
      min_days_per_season = min_days_per_season,
      a_bounds = a_bounds,
      n_warm_bounds = n_warm_bounds,
      n_cold_bounds = n_cold_bounds,
      tau_bounds = tau_bounds,
      initial_parameters = if (is_permafrost_group) {
        initial_parameters
      } else {
        NULL
      },
      maxit = maxit,
      output_dir = group_output_dir
    )
  }

  # Use `$` explicitly so the logical vector is resolved from lookup_dt rather
  # than searched for in the calling scope by recent data.table versions.
  non_permafrost_lookup <- lookup_dt[
    lookup_dt$is_permafrost == FALSE
  ]
  permafrost_lookup <- lookup_dt[
    lookup_dt$is_permafrost == TRUE
  ]

  message(
    "NEON SoilT MLE classification: ",
    nrow(non_permafrost_lookup),
    " non-permafrost site(s), ",
    nrow(permafrost_lookup),
    " permafrost site(s); latitude threshold = ",
    latitude_threshold,
    "."
  )

  non_permafrost_output_dir <- if (is.null(output_dir)) {
    NULL
  } else {
    file.path(output_dir, "non_permafrost")
  }
  permafrost_output_dir <- if (is.null(output_dir)) {
    NULL
  } else {
    file.path(output_dir, "permafrost")
  }

  non_permafrost_fit <- run_group(
    non_permafrost_lookup,
    "non_permafrost",
    non_permafrost_output_dir
  )
  permafrost_fit <- run_group(
    permafrost_lookup,
    "permafrost",
    permafrost_output_dir
  )

  add_classification <- function(table, group_lookup) {
    if (is.null(table) || nrow(table) == 0L) {
      return(data.table::data.table())
    }

    merge(
      data.table::as.data.table(data.table::copy(table)),
      group_lookup[
        ,
        .(index, site_lat, is_permafrost, process_model)
      ],
      by = "index",
      all.x = TRUE,
      sort = FALSE
    )
  }

  best_tau <- data.table::rbindlist(
    list(
      add_classification(
        if (is.null(non_permafrost_fit)) NULL else non_permafrost_fit$summary,
        non_permafrost_lookup
      ),
      add_classification(
        if (is.null(permafrost_fit)) NULL else permafrost_fit$summary,
        permafrost_lookup
      )
    ),
    use.names = TRUE,
    fill = TRUE
  )

  if (nrow(best_tau) > 0L) {
    data.table::setorder(best_tau, index, verticalPosition)
  }

  select_loyo <- function(fit) {
    if (is.null(fit)) {
      return(empty_loyo_result())
    }

    list(
      year_coverage = fit$year_coverage,
      fold_parameters = fit$fold_parameters,
      fold_metrics = fit$fold_metrics,
      validation = fit$validation,
      status = fit$status
    )
  }

  result <- list(
    best_tau = best_tau,
    non_permafrost_loyo = select_loyo(non_permafrost_fit),
    permafrost_loyo = select_loyo(permafrost_fit)
  )

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(
      best_tau,
      file.path(output_dir, "neon_tsoil_best_tau.csv")
    )
    saveRDS(
      result,
      file.path(output_dir, "neon_tsoil_mle_results.rds"),
      compress = FALSE
    )
  }

  return(result)
}


# Example:
# neon_mle_results <- neon_tsoil_mle(
#   lookup = lookup,
#   start_year = 2017L,
#   end_year = 2024L,
#   multi_site = multi_site,
#   era5_all = era5_all,
#   vertical_positions = "502",
#   workers = 5L,
#   permafrost_latitude_threshold = 63,
#   permafrost_fit_tau_only = TRUE,
#   permafrost_fixed_a_C = 0,
#   permafrost_fixed_n_warm = 1,
#   permafrost_fixed_n_cold = 0.3,
#   tau_bounds = c(0.125, 180),
#   warmup_days = 180L,
#   min_observations = 100L,
#   min_days_per_year = 120L,
#   min_days_per_season = 10L,
#   output_dir = "/path/to/neon_tsoil_mle"
# )
