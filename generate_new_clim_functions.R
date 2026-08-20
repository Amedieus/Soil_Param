#' Generate soil temperature with a causal exponential filter
#'
#' Applies the first-order soil-temperature response used by
#' [met2model.SIPNET_tau()]. For timestep \eqn{t},
#'
#' \deqn{\alpha_t = 1 - \exp(-\Delta t_t / \tau)}
#'
#' \deqn{T_{soil,t} = T_{soil,t-1} +
#'   \alpha_t (T_{air,t} - T_{soil,t-1}).}
#'
#' Unlike `stats::convolve()`, this recursion is causal and does not wrap the
#' end of a time series back to its beginning. A supplied initial state makes
#' it possible to carry soil temperature between consecutive files or years.
#'
#' @param Tair_C Numeric vector of air temperature in degrees Celsius.
#' @param dt_days Numeric timestep length in days. May be a scalar or have the
#'   same length as `Tair_C`.
#' @param tau_days Positive soil-temperature response time in days.
#' @param initial_soilT_C Optional soil-temperature state immediately before
#'   the first timestep. If `NULL`, the first air temperature is used.
#'
#' @return A list containing `soilT_C`, `final_soilT_C`, and the timestep-level
#'   response coefficient `alpha`.
#'
#' @md
#' @export
#' @author Yang Gu
generate_soilT_from_tau <- function(Tair_C,
                                    dt_days,
                                    tau_days,
                                    initial_soilT_C = NULL) {
  Tair_C <- as.numeric(Tair_C)
  dt_days <- as.numeric(dt_days)
  tau_days <- as.numeric(tau_days)
  
  if (length(Tair_C) == 0L) {
    return(list(
      soilT_C = numeric(),
      final_soilT_C = NA_real_,
      alpha = numeric()
    ))
  }
  if (any(!is.finite(Tair_C))) {
    stop("`Tair_C` must contain only finite values.", call. = FALSE)
  }
  if (length(dt_days) == 1L) {
    dt_days <- rep(dt_days, length(Tair_C))
  }
  if (length(dt_days) != length(Tair_C) ||
      any(!is.finite(dt_days)) ||
      any(dt_days <= 0)) {
    stop(
      "`dt_days` must be positive and have length 1 or length(Tair_C).",
      call. = FALSE
    )
  }
  if (length(tau_days) != 1L ||
      !is.finite(tau_days) ||
      tau_days <= 0) {
    stop("`tau_days` must be one positive finite value.", call. = FALSE)
  }
  
  if (is.null(initial_soilT_C)) {
    state <- Tair_C[[1L]]
  } else {
    initial_soilT_C <- as.numeric(initial_soilT_C)
    if (length(initial_soilT_C) != 1L || !is.finite(initial_soilT_C)) {
      stop(
        "`initial_soilT_C` must be NULL or one finite value.",
        call. = FALSE
      )
    }
    state <- initial_soilT_C
  }
  
  alpha <- -expm1(-dt_days / tau_days)
  soilT_C <- numeric(length(Tair_C))
  soilT_C[[1L]] <- state
  
  if (length(Tair_C) > 1L) {
    for (i in seq.int(2L, length(Tair_C))) {
      state <- state + alpha[[i]] * (Tair_C[[i]] - state)
      soilT_C[[i]] <- state
    }
  }
  
  list(
    soilT_C = soilT_C,
    final_soilT_C = soilT_C[[length(soilT_C)]],
    alpha = alpha
  )
}


#' Calculate SIPNET air and soil vapor pressure deficits
#'
#' Calculates atmospheric vapor pressure and the vapor pressure deficit used
#' by SIPNET for soil evaporation. The soil deficit is
#'
#' \deqn{VPD_{soil} = e_s(T_{soil}) - e_a,}
#'
#' where \eqn{e_s(T_{soil})} is saturation vapor pressure evaluated at the
#' modeled soil temperature and \eqn{e_a} is canopy-air vapor pressure.
#'
#' Atmospheric vapor pressure can be supplied directly, recovered from an
#' input air VPD, or calculated from specific humidity and air pressure using
#' the same PEcAn atmospheric utilities used by `met2model.SIPNET()`.
#'
#' @param Tair_C Numeric vector of air temperature in degrees Celsius.
#' @param soilT_C Numeric vector of soil temperature in degrees Celsius.
#' @param air_vapor_pressure_Pa Optional canopy-air vapor pressure in Pa.
#' @param specific_humidity Optional specific humidity in kg kg-1. Required
#'   only when neither `air_vapor_pressure_Pa` nor `VPDair_input_Pa` is given.
#' @param air_pressure_Pa Optional air pressure in Pa. Required with
#'   `specific_humidity`.
#' @param VPDair_input_Pa Optional atmospheric VPD in Pa.
#' @param clamp_soil_vpd Logical. If `TRUE`, negative soil VPD is set to zero,
#'   matching the nonnegative relative-humidity treatment in PEcAn.
#'
#' @return A list containing saturation vapor pressures, canopy-air vapor
#'   pressure, air VPD, and soil VPD, all in Pa.
#'
#' @md
#' @export
#' @author Yang Gu
calculate_sipnet_vpd <- function(Tair_C,
                                 soilT_C,
                                 air_vapor_pressure_Pa = NULL,
                                 specific_humidity = NULL,
                                 air_pressure_Pa = NULL,
                                 VPDair_input_Pa = NULL,
                                 clamp_soil_vpd = TRUE) {
  recycle_numeric <- function(x, n, name) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, n)
    }
    if (length(x) != n || any(!is.finite(x))) {
      stop(
        "`", name, "` must be finite and have length 1 or ", n, ".",
        call. = FALSE
      )
    }
    x
  }
  
  Tair_C <- as.numeric(Tair_C)
  soilT_C <- as.numeric(soilT_C)
  n <- length(Tair_C)
  
  if (n == 0L || length(soilT_C) != n ||
      any(!is.finite(Tair_C)) || any(!is.finite(soilT_C))) {
    stop(
      "`Tair_C` and `soilT_C` must be finite vectors of equal, nonzero length.",
      call. = FALSE
    )
  }
  if (!requireNamespace("PEcAn.data.atmosphere", quietly = TRUE) ||
      !requireNamespace("PEcAn.utils", quietly = TRUE)) {
    stop(
      "Packages `PEcAn.data.atmosphere` and `PEcAn.utils` are required.",
      call. = FALSE
    )
  }
  
  saturation_air_Pa <- PEcAn.utils::ud_convert(
    PEcAn.data.atmosphere::get.es(Tair_C),
    "millibar",
    "Pa"
  )
  saturation_soil_Pa <- PEcAn.utils::ud_convert(
    PEcAn.data.atmosphere::get.es(soilT_C),
    "millibar",
    "Pa"
  )
  
  if (!is.null(air_vapor_pressure_Pa)) {
    air_vapor_pressure_Pa <- recycle_numeric(
      air_vapor_pressure_Pa,
      n,
      "air_vapor_pressure_Pa"
    )
    VPDair_Pa <- saturation_air_Pa - air_vapor_pressure_Pa
  } else if (!is.null(VPDair_input_Pa)) {
    VPDair_Pa <- recycle_numeric(VPDair_input_Pa, n, "VPDair_input_Pa")
    air_vapor_pressure_Pa <- saturation_air_Pa - VPDair_Pa
  } else {
    if (is.null(specific_humidity) || is.null(air_pressure_Pa)) {
      stop(
        "Supply canopy-air vapor pressure, air VPD, or both specific humidity and air pressure.",
        call. = FALSE
      )
    }
    specific_humidity <- recycle_numeric(
      specific_humidity,
      n,
      "specific_humidity"
    )
    air_pressure_Pa <- recycle_numeric(
      air_pressure_Pa,
      n,
      "air_pressure_Pa"
    )
    relative_humidity <- PEcAn.data.atmosphere::qair2rh(
      specific_humidity,
      Tair_C,
      press = air_pressure_Pa / 100
    )
    VPDair_Pa <- saturation_air_Pa * (1 - relative_humidity)
    air_vapor_pressure_Pa <- saturation_air_Pa - VPDair_Pa
  }
  
  VPDsoil_Pa <- saturation_soil_Pa - air_vapor_pressure_Pa
  if (isTRUE(clamp_soil_vpd)) {
    VPDsoil_Pa <- pmax(VPDsoil_Pa, 0)
  }
  
  list(
    saturation_air_Pa = saturation_air_Pa,
    saturation_soil_Pa = saturation_soil_Pa,
    air_vapor_pressure_Pa = air_vapor_pressure_Pa,
    VPDair_Pa = VPDair_Pa,
    VPDsoil_Pa = VPDsoil_Pa
  )
}

.sipnet_tau_log <- function(level, ...) {
  if (requireNamespace("PEcAn.logger", quietly = TRUE)) {
    logger <- getExportedValue("PEcAn.logger", paste0("logger.", level))
    do.call(logger, list(...))
  } else {
    message(paste0(..., collapse = ""))
  }
}

.sipnet_tau_filename <- function(in.prefix,
                                 start_date,
                                 end_date,
                                 year.fragment) {
  if (isTRUE(year.fragment)) {
    if (grepl("\\.nc$", in.prefix)) {
      return(sub("\\.nc$", ".clim", in.prefix))
    }
    return(paste0(in.prefix, ".clim"))
  }
  
  paste(
    in.prefix,
    format(as.Date(start_date), "%Y-%m-%d"),
    format(as.Date(end_date), "%Y-%m-%d"),
    "clim",
    sep = "."
  )
}

.sipnet_tau_result <- function(file,
                               start_date,
                               end_date,
                               soil_tau_days,
                               tau_applied,
                               status) {
  host <- if (requireNamespace("PEcAn.remote", quietly = TRUE)) {
    PEcAn.remote::fqdn()
  } else {
    unname(Sys.info()[["nodename"]])
  }
  
  data.frame(
    file = file,
    host = host,
    mimetype = "text/csv",
    formatname = "Sipnet.climna",
    startdate = as.POSIXlt(start_date, tz = "UTC"),
    enddate = as.POSIXlt(end_date, tz = "UTC"),
    dbfile.name = basename(file),
    soil_tau_days = soil_tau_days,
    soil_temperature_source = "tau",
    tau_applied = tau_applied,
    status = status,
    stringsAsFactors = FALSE
  )
}


.write_sipnet_clim_atomic <- function(clim, path, overwrite) {
  output_dir <- dirname(path)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = output_dir
  )
  on.exit(unlink(temporary_file, force = TRUE), add = TRUE)
  
  utils::write.table(
    x = format(as.data.frame(clim), digits = 4),
    file = temporary_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Output exists and `overwrite = FALSE`: ", path, call. = FALSE)
  }
  
  backup_file <- NULL
  if (file.exists(path)) {
    backup_file <- tempfile(
      pattern = paste0(basename(path), ".backup."),
      tmpdir = output_dir
    )
    if (!file.rename(path, backup_file)) {
      stop("Could not move the existing output before replacement: ", path)
    }
  }
  
  installed <- file.rename(temporary_file, path)
  if (!installed) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }
    stop("Could not install the completed SIPNET climate file: ", path)
  }
  
  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file, force = TRUE)
  }
  
  invisible(path)
}


#' Generate a tau-corrected SIPNET climate file
#'
#' Uses the official [PEcAn.SIPNET::met2model.SIPNET()] function to generate a
#' standard SIPNET climate file, then replaces soil temperature with a causal
#' tau-filtered estimate and regenerates soil VPD from the new soil temperature
#' and PEcAn-generated canopy-air vapor pressure. All other forcing columns are
#' retained exactly as produced by PEcAn.
#'
#' SIPNET v2 climate files contain 12 columns. Legacy v1 files contain the same
#' 12 forcing columns plus a leading grid index and trailing soil-wetness
#' placeholder. This function detects both layouts and updates the appropriate
#' soil-temperature and soil-VPD columns. It never changes only soil
#' temperature while leaving the dependent soil VPD unchanged.
#'
#' The base PEcAn file is generated in a temporary directory. The completed
#' tau-corrected file is validated for missing and non-finite values before it
#' is installed in `outfolder`.
#'
#' @inheritParams PEcAn.SIPNET::met2model.SIPNET
#' @param soil_tau_days Positive soil-temperature response time in days.
#' @param initial_soilT_C Optional soil-temperature state immediately before
#'   the first output timestep. If `NULL`, the first air temperature is used.
#' @param clamp_soil_vpd Logical. If `TRUE`, negative soil VPD is set to zero.
#' @param met2model_function Optional function used to generate the base SIPNET
#'   climate file. `NULL` uses `PEcAn.SIPNET::met2model.SIPNET()`. This argument
#'   is primarily intended for testing.
#'
#' @return Invisibly returns a data.frame describing the output file, tau, and
#'   generation status.
#'
#' @references
#' PEcAn SIPNET meteorological converter:
#' <https://github.com/PecanProject/pecan/blob/develop/models/sipnet/R/met2model.SIPNET.R>
#'
#' SIPNET climate-input documentation:
#' <https://github.com/PecanProject/sipnet/blob/master/docs/user-guide/model-inputs.md>
#'
#' @md
#' @export
#' @author Yang Gu
met2model.SIPNET_tau <- function(in.path,
                                 in.prefix,
                                 outfolder,
                                 start_date,
                                 end_date,
                                 soil_tau_days,
                                 var.names = NULL,
                                 initial_soilT_C = NULL,
                                 overwrite = FALSE,
                                 verbose = FALSE,
                                 year.fragment = FALSE,
                                 clim_format_version = c("v2", "v1"),
                                 clamp_soil_vpd = TRUE,
                                 met2model_function = NULL,
                                 ...) {
  clim_format_version <- match.arg(clim_format_version)
  soil_tau_days <- as.numeric(soil_tau_days)
  
  if (length(soil_tau_days) != 1L ||
      !is.finite(soil_tau_days) ||
      soil_tau_days <= 0) {
    stop("`soil_tau_days` must be one positive finite value.", call. = FALSE)
  }
  start_date_parsed <- as.Date(start_date)
  end_date_parsed <- as.Date(end_date)
  if (length(start_date_parsed) != 1L ||
      length(end_date_parsed) != 1L ||
      is.na(start_date_parsed) ||
      is.na(end_date_parsed) ||
      end_date_parsed < start_date_parsed) {
    stop(
      "`start_date` and `end_date` must be valid dates with end >= start.",
      call. = FALSE
    )
  }
  
  output_file <- file.path(
    outfolder,
    .sipnet_tau_filename(
      in.prefix,
      start_date_parsed,
      end_date_parsed,
      year.fragment
    )
  )
  
  if (file.exists(output_file) && !isTRUE(overwrite)) {
    if (isTRUE(verbose)) {
      .sipnet_tau_log(
        "info",
        "File already exists; tau forcing was not regenerated: ",
        output_file
      )
    }
    return(invisible(.sipnet_tau_result(
      file = output_file,
      start_date = start_date_parsed,
      end_date = end_date_parsed,
      soil_tau_days = soil_tau_days,
      tau_applied = FALSE,
      status = "EXISTS_NOT_OVERWRITTEN"
    )))
  }
  
  if (is.null(met2model_function)) {
    if (!requireNamespace("PEcAn.SIPNET", quietly = TRUE)) {
      stop(
        "Package `PEcAn.SIPNET` is required unless `met2model_function` is supplied.",
        call. = FALSE
      )
    }
    met2model_function <- getExportedValue(
      "PEcAn.SIPNET",
      "met2model.SIPNET"
    )
  }
  if (!is.function(met2model_function)) {
    stop("`met2model_function` must be NULL or a function.", call. = FALSE)
  }
  
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  temporary_outfolder <- tempfile(
    pattern = "sipnet_tau_base_",
    tmpdir = outfolder
  )
  dir.create(temporary_outfolder, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temporary_outfolder, recursive = TRUE, force = TRUE), add = TRUE)
  
  if (isTRUE(verbose)) {
    .sipnet_tau_log(
      "info",
      "Generating PEcAn base climate forcing before applying tau = ",
      soil_tau_days,
      " days"
    )
  }
  
  base_result <- met2model_function(
    in.path = in.path,
    in.prefix = in.prefix,
    outfolder = temporary_outfolder,
    start_date = start_date,
    end_date = end_date,
    var.names = var.names,
    overwrite = TRUE,
    verbose = verbose,
    year.fragment = year.fragment,
    clim_format_version = clim_format_version,
    ...
  )
  
  if (is.null(base_result) ||
      !"file" %in% names(base_result) ||
      length(base_result$file) == 0L ||
      !file.exists(base_result$file[[1L]])) {
    stop("PEcAn did not produce a base SIPNET climate file.", call. = FALSE)
  }
  
  base_file <- base_result$file[[1L]]
  output_file <- file.path(outfolder, basename(base_file))
  clim <- utils::read.table(
    base_file,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (nrow(clim) == 0L) {
    stop("The PEcAn base climate file contains no rows.", call. = FALSE)
  }
  if (ncol(clim) == 12L) {
    timestep_column <- 4L
    air_temperature_column <- 5L
    soil_temperature_column <- 6L
    soil_vpd_column <- 10L
    air_vapor_pressure_column <- 11L
  } else if (ncol(clim) == 14L) {
    timestep_column <- 5L
    air_temperature_column <- 6L
    soil_temperature_column <- 7L
    soil_vpd_column <- 11L
    air_vapor_pressure_column <- 12L
  } else {
    stop(
      "Unrecognized SIPNET climate format: expected 12 or 14 columns, found ",
      ncol(clim),
      ".",
      call. = FALSE
    )
  }
  
  clim[] <- lapply(clim, function(x) suppressWarnings(as.numeric(x)))
  if (anyNA(clim)) {
    stop(
      "The PEcAn base climate file contains nonnumeric or missing values.",
      call. = FALSE
    )
  }
  
  timestep_days <- clim[[timestep_column]]
  timestep_days[timestep_days < 0] <-
    -timestep_days[timestep_days < 0] / 86400
  
  soil_temperature <- generate_soilT_from_tau(
    Tair_C = clim[[air_temperature_column]],
    dt_days = timestep_days,
    tau_days = soil_tau_days,
    initial_soilT_C = initial_soilT_C
  )
  vapor_pressure <- calculate_sipnet_vpd(
    Tair_C = clim[[air_temperature_column]],
    soilT_C = soil_temperature$soilT_C,
    air_vapor_pressure_Pa = clim[[air_vapor_pressure_column]],
    clamp_soil_vpd = clamp_soil_vpd
  )
  
  clim[[soil_temperature_column]] <- soil_temperature$soilT_C
  clim[[soil_vpd_column]] <- vapor_pressure$VPDsoil_Pa
  
  clim_matrix <- as.matrix(clim)
  if (anyNA(clim_matrix) || any(!is.finite(clim_matrix))) {
    stop(
      "Tau correction produced missing or non-finite SIPNET forcing values.",
      call. = FALSE
    )
  }
  
  .write_sipnet_clim_atomic(
    clim = clim,
    path = output_file,
    overwrite = overwrite
  )
  
  result <- base_result
  result$file <- output_file
  result$dbfile.name <- basename(output_file)
  result$soil_tau_days <- soil_tau_days
  result$soil_temperature_source <- "tau"
  result$initial_soilT_C <- if (is.null(initial_soilT_C)) {
    NA_real_
  } else {
    as.numeric(initial_soilT_C)
  }
  result$final_soilT_C <- soil_temperature$final_soilT_C
  result$tau_applied <- TRUE
  result$status <- "GENERATED"
  
  if (isTRUE(verbose)) {
    .sipnet_tau_log("info", "Wrote tau-corrected forcing: ", output_file)
  }
  
  invisible(result)
}


.extract_pft_tau_table <- function(pft_tau_test) {
  if (is.data.frame(pft_tau_test)) {
    pft_tau <- pft_tau_test
  } else if (is.list(pft_tau_test) &&
             !is.null(pft_tau_test$pft_tau) &&
             is.data.frame(pft_tau_test$pft_tau)) {
    pft_tau <- pft_tau_test$pft_tau
  } else {
    stop(
      "`pft_tau_test` must be a data.frame or contain data.frame `pft_tau`.",
      call. = FALSE
    )
  }
  
  required_columns <- c("final_pft", "pft_tau_days")
  missing_columns <- setdiff(required_columns, names(pft_tau))
  if (length(missing_columns) > 0L) {
    stop(
      "PFT tau table is missing: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  pft_tau
}


.safe_pft_name <- function(x) {
  safe <- gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(x)))
  safe <- gsub("^_+|_+$", "", safe)
  if (!nzchar(safe)) {
    stop("PFT name cannot be converted to a safe file name.", call. = FALSE)
  }
  safe
}


#' Generate tau-corrected SIPNET climate files for one PFT
#'
#' Selects one PFT-level tau, finds every model index assigned to that PFT,
#' and generates SIPNET forcing for each requested ERA5 ensemble member. Input
#' directories are expected to follow
#' `ERA5_<index>_<member>/ERA5.<member>.<year>.nc`. Output directories preserve
#' the `ERA5_<index>_<member>` layout so that they can replace the meteorology
#' root in a subsequent PEcAn/SIPNET run.
#'
#' Jobs are parallelized over index-member combinations on Unix-like systems.
#' Every attempted job is recorded in a PFT-specific manifest, including
#' missing inputs and generation errors.
#'
#' @param lookup A data.frame or data.table containing unique model `index`
#'   values.
#' @param newpft A data.frame or data.table containing `index` and `final_pft`.
#' @param pft_tau_test A data.frame with `final_pft` and `pft_tau_days`, or a
#'   list containing that table as `pft_tau`.
#' @param pft_name Character name of the PFT to process.
#' @param input_root Root containing PEcAn-standardized ERA5 directories.
#' @param output_root Root where tau-corrected forcing directories and the
#'   manifest will be written.
#' @param members Integer ERA5 ensemble members to process.
#' @param start_date First forcing date, inclusive, in `YYYY-MM-DD` format.
#' @param end_date Last forcing date, inclusive, in `YYYY-MM-DD` format.
#' @param n_cores Number of forked workers on Unix-like systems.
#' @param overwrite Logical indicating whether existing tau climate files may
#'   be replaced.
#' @param verbose Logical indicating whether to print progress messages.
#' @param clim_format_version SIPNET climate format. Current SIPNET uses `"v2"`
#'   (12 columns); `"v1"` writes the legacy 14-column format.
#' @param initial_soilT_C Optional initial soil-temperature state passed to
#'   [met2model.SIPNET_tau()].
#' @param clamp_soil_vpd Logical indicating whether negative soil VPD is set
#'   to zero.
#' @param stop_on_error Logical. If `TRUE`, stop after writing the manifest when
#'   any job failed or had missing input.
#' @param met2model_function Optional base meteorological converter passed to
#'   [met2model.SIPNET_tau()].
#'
#' @return Invisibly returns the job manifest as a data.frame.
#'
#' @md
#' @export
#' @author Yang Gu
generate_tau_clims_for_pft <- function(
    lookup,
    newpft,
    pft_tau_test,
    pft_name,
    input_root =
      "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024",
    output_root =
      "/projectnb/dietzelab/guYANG/pecan/modified_met/pft_specific_tau",
    members = 1:10,
    start_date = "2012-01-01",
    end_date = "2024-12-31",
    n_cores = 6L,
    overwrite = TRUE,
    verbose = TRUE,
    clim_format_version = c("v2", "v1"),
    initial_soilT_C = NULL,
    clamp_soil_vpd = TRUE,
    stop_on_error = FALSE,
    met2model_function = NULL) {
  clim_format_version <- match.arg(clim_format_version)
  lookup <- as.data.frame(lookup, stringsAsFactors = FALSE)
  newpft <- as.data.frame(newpft, stringsAsFactors = FALSE)
  pft_tau <- .extract_pft_tau_table(pft_tau_test)
  
  if (!"index" %in% names(lookup)) {
    stop("`lookup` must contain `index`.", call. = FALSE)
  }
  if (!all(c("index", "final_pft") %in% names(newpft))) {
    stop("`newpft` must contain `index` and `final_pft`.", call. = FALSE)
  }
  if (anyDuplicated(lookup$index)) {
    stop("`lookup$index` must be unique.", call. = FALSE)
  }
  if (anyDuplicated(newpft$index)) {
    stop("`newpft$index` must be unique.", call. = FALSE)
  }
  
  pft_name <- as.character(pft_name)
  if (length(pft_name) != 1L || is.na(pft_name) || !nzchar(pft_name)) {
    stop("`pft_name` must be one nonempty character value.", call. = FALSE)
  }
  
  tau_rows <- pft_tau[as.character(pft_tau$final_pft) == pft_name, , drop = FALSE]
  if (nrow(tau_rows) != 1L) {
    stop(
      "Expected exactly one PFT tau row for `", pft_name,
      "`; found ", nrow(tau_rows), ".",
      call. = FALSE
    )
  }
  tau_days <- as.numeric(tau_rows$pft_tau_days[[1L]])
  if (!is.finite(tau_days) || tau_days <= 0) {
    stop("PFT tau must be one positive finite value.", call. = FALSE)
  }
  
  pft_indices <- unique(newpft$index[
    as.character(newpft$final_pft) == pft_name
  ])
  pft_indices <- suppressWarnings(as.integer(pft_indices))
  pft_indices <- pft_indices[is.finite(pft_indices)]
  if (length(pft_indices) == 0L) {
    stop("No model indices are assigned to PFT `", pft_name, "`.", call. = FALSE)
  }
  
  missing_lookup_indices <- setdiff(pft_indices, as.integer(lookup$index))
  if (length(missing_lookup_indices) > 0L) {
    stop(
      "PFT indices missing from `lookup`: ",
      paste(missing_lookup_indices, collapse = ", "),
      call. = FALSE
    )
  }
  
  members <- unique(suppressWarnings(as.integer(members)))
  if (length(members) == 0L || any(!is.finite(members))) {
    stop("`members` must contain valid integers.", call. = FALSE)
  }
  n_cores <- as.integer(n_cores)
  if (length(n_cores) != 1L || !is.finite(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be one positive integer.", call. = FALSE)
  }
  
  start_year <- as.integer(format(as.Date(start_date), "%Y"))
  end_year <- as.integer(format(as.Date(end_date), "%Y"))
  if (!is.finite(start_year) || !is.finite(end_year) || end_year < start_year) {
    stop("Invalid forcing date range.", call. = FALSE)
  }
  years <- seq.int(start_year, end_year)
  jobs <- expand.grid(
    index = sort(pft_indices),
    member = sort(members),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  safe_pft <- .safe_pft_name(pft_name)
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  
  run_job <- function(i) {
    index_i <- jobs$index[[i]]
    member_i <- jobs$member[[i]]
    input_dir <- file.path(
      input_root,
      sprintf("ERA5_%d_%d", index_i, member_i)
    )
    output_dir <- file.path(
      output_root,
      sprintf("ERA5_%d_%d", index_i, member_i)
    )
    met_prefix <- sprintf("ERA5.%d", member_i)
    expected_inputs <- file.path(
      input_dir,
      sprintf("%s.%d.nc", met_prefix, years)
    )
    missing_inputs <- expected_inputs[!file.exists(expected_inputs)]
    
    base_manifest <- data.frame(
      index = index_i,
      member = member_i,
      final_pft = pft_name,
      tau_days = tau_days,
      input_dir = input_dir,
      output_dir = output_dir,
      new_clim_path = NA_character_,
      status = NA_character_,
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (length(missing_inputs) > 0L) {
      base_manifest$status <- "MISSING_INPUT"
      base_manifest$error_message <- paste(
        "Missing",
        length(missing_inputs),
        "annual NetCDF file(s); first:",
        missing_inputs[[1L]]
      )
      return(base_manifest)
    }
    
    tryCatch({
      result <- met2model.SIPNET_tau(
        in.path = input_dir,
        in.prefix = met_prefix,
        outfolder = output_dir,
        start_date = start_date,
        end_date = end_date,
        soil_tau_days = tau_days,
        initial_soilT_C = initial_soilT_C,
        overwrite = overwrite,
        verbose = FALSE,
        clim_format_version = clim_format_version,
        clamp_soil_vpd = clamp_soil_vpd,
        met2model_function = met2model_function
      )
      
      if (is.null(result) || !"file" %in% names(result)) {
        stop("No climate-file result was returned.")
      }
      base_manifest$new_clim_path <- as.character(result$file[[1L]])
      base_manifest$status <- if ("status" %in% names(result)) {
        as.character(result$status[[1L]])
      } else {
        "GENERATED"
      }
      base_manifest
    }, error = function(e) {
      base_manifest$status <- "ERROR"
      base_manifest$error_message <- conditionMessage(e)
      base_manifest
    })
  }
  
  if (isTRUE(verbose)) {
    .sipnet_tau_log(
      "info",
      "Generating ",
      nrow(jobs),
      " tau climate file(s) for PFT `",
      pft_name,
      "` with tau = ",
      tau_days,
      " days"
    )
  }
  
  job_ids <- seq_len(nrow(jobs))
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    job_results <- parallel::mclapply(
      job_ids,
      run_job,
      mc.cores = min(n_cores, length(job_ids)),
      mc.preschedule = FALSE
    )
  } else {
    if (n_cores > 1L && .Platform$OS.type == "windows" && isTRUE(verbose)) {
      .sipnet_tau_log(
        "warn",
        "Forked processing is unavailable on Windows; running sequentially."
      )
    }
    job_results <- lapply(job_ids, run_job)
  }
  
  manifest <- do.call(rbind, job_results)
  rownames(manifest) <- NULL
  manifest <- manifest[order(manifest$index, manifest$member), , drop = FALSE]
  manifest_path <- file.path(
    output_root,
    paste0(safe_pft, "_clim_manifest.csv")
  )
  utils::write.csv(manifest, manifest_path, row.names = FALSE, na = "")
  
  if (isTRUE(verbose)) {
    .sipnet_tau_log("info", "Wrote climate manifest: ", manifest_path)
  }
  
  failed <- manifest$status %in% c("MISSING_INPUT", "ERROR")
  if (isTRUE(stop_on_error) && any(failed)) {
    stop(
      sum(failed),
      " climate job(s) failed. See ",
      manifest_path,
      ".",
      call. = FALSE
    )
  }
  
  invisible(manifest)
}


#' Generate tau-corrected SIPNET climate files for all PFTs
#'
#' Calls [generate_tau_clims_for_pft()] for each requested PFT. Processing is
#' sequential across PFTs because the single-PFT function already parallelizes
#' index-member jobs. PFT-level manifests and one combined manifest are written
#' under `output_root`.
#'
#' @inheritParams generate_tau_clims_for_pft
#' @param pfts Optional character vector of PFT names. `NULL` processes every
#'   PFT represented in both `newpft` and the PFT tau table.
#' @param missing_tau_action Action when a requested PFT has no unique finite
#'   positive tau: `"error"` or `"skip"`.
#'
#' @return Invisibly returns the combined job manifest as a data.frame.
#'
#' @md
#' @export
#' @author Yang Gu
generate_tau_clims_for_all_pfts <- function(
    lookup,
    newpft,
    pft_tau_test,
    pfts = NULL,
    input_root =
      "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024",
    output_root =
      "/projectnb/dietzelab/guYANG/pecan/modified_met/pft_tau_all",
    members = 1:10,
    start_date = "2012-01-01",
    end_date = "2024-12-31",
    n_cores = 18L,
    overwrite = TRUE,
    verbose = TRUE,
    missing_tau_action = c("error", "skip"),
    clim_format_version = c("v2", "v1"),
    initial_soilT_C = NULL,
    clamp_soil_vpd = TRUE,
    stop_on_error = FALSE,
    met2model_function = NULL) {
  missing_tau_action <- match.arg(missing_tau_action)
  clim_format_version <- match.arg(clim_format_version)
  newpft <- as.data.frame(newpft, stringsAsFactors = FALSE)
  pft_tau <- .extract_pft_tau_table(pft_tau_test)
  
  if (!all(c("index", "final_pft") %in% names(newpft))) {
    stop("`newpft` must contain `index` and `final_pft`.", call. = FALSE)
  }
  
  assigned_pfts <- unique(as.character(newpft$final_pft))
  assigned_pfts <- assigned_pfts[!is.na(assigned_pfts) & nzchar(assigned_pfts)]
  tau_pfts <- unique(as.character(pft_tau$final_pft))
  
  if (is.null(pfts)) {
    pfts <- sort(intersect(assigned_pfts, tau_pfts))
  } else {
    pfts <- unique(as.character(pfts))
    pfts <- pfts[!is.na(pfts) & nzchar(pfts)]
  }
  if (length(pfts) == 0L) {
    stop("No PFTs are available for climate generation.", call. = FALSE)
  }
  
  tau_count <- vapply(pfts, function(pft_name) {
    rows <- pft_tau[as.character(pft_tau$final_pft) == pft_name, , drop = FALSE]
    tau_value <- suppressWarnings(as.numeric(rows$pft_tau_days))
    as.integer(
      nrow(rows) == 1L &&
        length(tau_value) == 1L &&
        is.finite(tau_value) &&
        tau_value > 0
    )
  }, integer(1))
  invalid_pfts <- pfts[tau_count != 1L]
  if (length(invalid_pfts) > 0L && missing_tau_action == "error") {
    stop(
      "PFTs without exactly one valid tau: ",
      paste(invalid_pfts, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(invalid_pfts) > 0L) {
    if (isTRUE(verbose)) {
      .sipnet_tau_log(
        "warn",
        "Skipping PFTs without exactly one valid tau: ",
        paste(invalid_pfts, collapse = ", ")
      )
    }
    pfts <- setdiff(pfts, invalid_pfts)
  }
  if (length(pfts) == 0L) {
    stop("No PFTs remain after tau validation.", call. = FALSE)
  }
  
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  manifests <- vector("list", length(pfts))
  
  for (i in seq_along(pfts)) {
    pft_name <- pfts[[i]]
    if (isTRUE(verbose)) {
      .sipnet_tau_log(
        "info",
        "[",
        i,
        "/",
        length(pfts),
        "] PFT `",
        pft_name,
        "`"
      )
    }
    
    manifests[[i]] <- tryCatch(
      generate_tau_clims_for_pft(
        lookup = lookup,
        newpft = newpft,
        pft_tau_test = pft_tau,
        pft_name = pft_name,
        input_root = input_root,
        output_root = output_root,
        members = members,
        start_date = start_date,
        end_date = end_date,
        n_cores = n_cores,
        overwrite = overwrite,
        verbose = verbose,
        clim_format_version = clim_format_version,
        initial_soilT_C = initial_soilT_C,
        clamp_soil_vpd = clamp_soil_vpd,
        stop_on_error = FALSE,
        met2model_function = met2model_function
      ),
      error = function(e) {
        data.frame(
          index = NA_integer_,
          member = NA_integer_,
          final_pft = pft_name,
          tau_days = NA_real_,
          input_dir = NA_character_,
          output_dir = output_root,
          new_clim_path = NA_character_,
          status = "PFT_ERROR",
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }
  
  combined_manifest <- do.call(rbind, manifests)
  rownames(combined_manifest) <- NULL
  combined_manifest_path <- file.path(
    output_root,
    "all_pfts_clim_manifest.csv"
  )
  utils::write.csv(
    combined_manifest,
    combined_manifest_path,
    row.names = FALSE,
    na = ""
  )
  
  if (isTRUE(verbose)) {
    .sipnet_tau_log(
      "info",
      "Wrote combined climate manifest: ",
      combined_manifest_path
    )
  }
  
  failed <- combined_manifest$status %in%
    c("MISSING_INPUT", "ERROR", "PFT_ERROR")
  if (isTRUE(stop_on_error) && any(failed)) {
    stop(
      sum(failed),
      " climate job(s) failed. See ",
      combined_manifest_path,
      ".",
      call. = FALSE
    )
  }
  
  invisible(combined_manifest)
}
