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
        "`SoilT_tau_permafrost_functions.R` before calling ",
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

  non_permafrost_lookup <- lookup_dt[!is_permafrost]
  permafrost_lookup <- lookup_dt[is_permafrost]

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
