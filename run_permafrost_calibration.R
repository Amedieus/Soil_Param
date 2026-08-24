#!/usr/bin/env Rscript

# =============================================================================
# Run the full permafrost soil-temperature calibration workflow
# =============================================================================
#
# This runner calls `run_all_permafrost_indices_depths_parallel()` from
# `SoilT_tau_functions.R`. For each requested site and NEON depth, the workflow
# estimates a_C, n_warm, n_cold, and tau_days, performs leave-one-year-out
# (LOYO) validation, and writes both site-level and combined output files.
#
# Command-line use:
#
# Rscript run_permafrost_calibration.R \
#   /path/to/permafrost_lookup.csv \
#   /projectnb/dietzelab/guYANG/soilparam/permafrost_tsoil \
#   16 \
#   501,502,503,504
#
# Interactive use:
#
# source("run_permafrost_calibration.R")
# permafrost_tsoil_results <- run_permafrost_calibration(
#   permafrost_lookup = permafrost,
#   output_root =
#     "/projectnb/dietzelab/guYANG/soilparam/permafrost_tsoil",
#   workers = 16L,
#   vertical_positions = c("501", "502", "503", "504")
# )
#
# The lookup must contain exactly one row per index and the columns `index` and
# `NEON_code`. `AmeriFlux_ID` is retained in the output when it is available.
#
# @author Yang Gu


# =============================================================================
# Locate this script and its companion function file
# =============================================================================

get_permafrost_runner_dir <- function() {
  command_args <- commandArgs(
    trailingOnly = FALSE
  )
  
  file_arg <- grep(
    "^--file=",
    command_args,
    value = TRUE
  )
  
  if (length(file_arg) == 0L) {
    return(
      normalizePath(
        getwd(),
        mustWork = TRUE
      )
    )
  }
  
  dirname(
    normalizePath(
      sub(
        "^--file=",
        "",
        file_arg[1L]
      ),
      mustWork = TRUE
    )
  )
}


PERMAFROST_RUNNER_DIR <- get_permafrost_runner_dir()


# =============================================================================
# Read a permafrost lookup supplied to the command-line runner
# =============================================================================

#' Read a permafrost-site lookup
#'
#' @param lookup_file Path to a CSV, TSV, RDS, RData, or rda file. An RData/rda
#'   file should contain an object named `permafrost`, `permafrost_lookup`, or
#'   `lookup`. If none of those names exists, the file must contain exactly one
#'   data-frame-like object.
#'
#' @return A `data.table` containing the permafrost lookup.
#'
#' @keywords internal
#' @author Yang Gu
read_permafrost_lookup <- function(
    lookup_file
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "Package `data.table` is required.",
      call. = FALSE
    )
  }
  
  lookup_file <- normalizePath(
    lookup_file,
    mustWork = TRUE
  )
  
  extension <- tolower(
    tools::file_ext(
      lookup_file
    )
  )
  
  if (extension == "csv") {
    lookup <- data.table::fread(
      lookup_file
    )
  } else if (extension %in% c("tsv", "txt")) {
    lookup <- data.table::fread(
      lookup_file,
      sep = "\t"
    )
  } else if (extension == "rds") {
    lookup <- readRDS(
      lookup_file
    )
  } else if (extension %in% c("rdata", "rda")) {
    lookup_environment <- new.env(
      parent = emptyenv()
    )
    
    loaded_names <- load(
      lookup_file,
      envir = lookup_environment
    )
    
    preferred_names <- c(
      "permafrost",
      "permafrost_lookup",
      "lookup"
    )
    
    selected_name <- intersect(
      preferred_names,
      loaded_names
    )
    
    if (length(selected_name) > 0L) {
      lookup <- get(
        selected_name[1L],
        envir = lookup_environment,
        inherits = FALSE
      )
    } else {
      data_frame_names <- loaded_names[
        vapply(
          loaded_names,
          function(object_name) {
            object <- get(
              object_name,
              envir = lookup_environment,
              inherits = FALSE
            )
            
            is.data.frame(object) ||
              data.table::is.data.table(object)
          },
          logical(1L)
        )
      ]
      
      if (length(data_frame_names) != 1L) {
        stop(
          paste0(
            "The RData/rda file must contain a named `permafrost` lookup ",
            "or exactly one data-frame-like object."
          ),
          call. = FALSE
        )
      }
      
      lookup <- get(
        data_frame_names[1L],
        envir = lookup_environment,
        inherits = FALSE
      )
    }
  } else {
    stop(
      paste0(
        "Unsupported lookup extension: .",
        extension,
        ". Use CSV, TSV, RDS, RData, or rda."
      ),
      call. = FALSE
    )
  }
  
  data.table::as.data.table(
    data.table::copy(
      lookup
    )
  )
}


# =============================================================================
# Validate the permafrost lookup before launching parallel jobs
# =============================================================================

#' Validate a permafrost-site lookup
#'
#' @param permafrost_lookup Data-frame-like lookup containing one row per model
#'   index and the columns `index` and `NEON_code`.
#'
#' @return A validated `data.table`.
#'
#' @keywords internal
#' @author Yang Gu
validate_permafrost_lookup <- function(
    permafrost_lookup
) {
  lookup_dt <- data.table::as.data.table(
    data.table::copy(
      permafrost_lookup
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
  
  if (length(missing_columns) > 0L) {
    stop(
      "Permafrost lookup is missing: ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  lookup_dt[
    ,
    index := suppressWarnings(
      as.integer(index)
    )
  ]
  
  lookup_dt[
    ,
    NEON_code := trimws(
      as.character(NEON_code)
    )
  ]
  
  invalid_rows <- lookup_dt[
    is.na(index) |
      is.na(NEON_code) |
      NEON_code == ""
  ]
  
  if (nrow(invalid_rows) > 0L) {
    stop(
      "Permafrost lookup contains missing/invalid `index` or `NEON_code`.",
      call. = FALSE
    )
  }
  
  duplicate_indices <- lookup_dt[
    ,
    .N,
    by = index
  ][
    N != 1L
  ]
  
  if (nrow(duplicate_indices) > 0L) {
    stop(
      paste0(
        "The permafrost workflow requires exactly one lookup row per index. ",
        "Duplicated indices: ",
        paste(
          duplicate_indices$index,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  
  data.table::setorder(
    lookup_dt,
    index
  )
  
  lookup_dt
}


# =============================================================================
# Public runner
# =============================================================================

#' Run full permafrost SoilT calibration and write result files
#'
#' @param permafrost_lookup Permafrost-only lookup with one row per `index`.
#' @param output_root Directory in which site-depth and combined outputs are
#'   written.
#' @param workers Number of parallel R workers.
#' @param vertical_positions NEON soil-temperature positions to process.
#' @param functions_file Path to `SoilT_tau_functions.R`.
#' @param obs_start_year First NEON observation year.
#' @param obs_end_year Last NEON observation year.
#' @param clim_root Root containing `ERA5_<index>_1` climate directories.
#' @param clim_basename SIPNET climate file name inside each climate directory.
#' @param neon_env_dir Root containing the NEON environmental data files.
#' @param ... Additional named arguments forwarded to
#'   `run_all_permafrost_indices_depths_parallel()` and then to
#'   `estimate_permafrost_tsoil_by_index()`.
#'
#' @return Invisibly, a list containing `summary` and `status`. The same object
#'   is also written to `permafrost_tsoil_results.rds` in `output_root`.
#'
#' @export
#' @author Yang Gu
run_permafrost_calibration <- function(
    permafrost_lookup,
    output_root =
      "/projectnb/dietzelab/guYANG/soilparam/permafrost_tsoil",
    workers = 16L,
    vertical_positions = c(
      "501",
      "502",
      "503",
      "504"
    ),
    functions_file = file.path(
      PERMAFROST_RUNNER_DIR,
      "SoilT_tau_functions.R"
    ),
    obs_start_year = 2017L,
    obs_end_year = 2024L,
    clim_root =
      "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024",
    clim_basename =
      "ERA5.1.2012-01-01.2024-12-31.clim",
    neon_env_dir =
      "/projectnb/dietzelab/jzobitz/02-NEON-sites/env-data",
    ...
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "Package `data.table` is required.",
      call. = FALSE
    )
  }
  
  if (!file.exists(functions_file)) {
    stop(
      "Cannot find permafrost function file: ",
      functions_file,
      call. = FALSE
    )
  }
  
  source(
    functions_file,
    local = globalenv()
  )
  
  if (!exists(
    "run_all_permafrost_indices_depths_parallel",
    mode = "function",
    inherits = TRUE
  )) {
    stop(
      paste0(
        "`run_all_permafrost_indices_depths_parallel()` was not loaded from ",
        functions_file,
        "."
      ),
      call. = FALSE
    )
  }
  
  lookup_dt <- validate_permafrost_lookup(
    permafrost_lookup
  )
  
  workers <- suppressWarnings(
    as.integer(workers)[1L]
  )
  
  if (is.na(workers) || workers < 1L) {
    stop(
      "`workers` must be one integer greater than or equal to one.",
      call. = FALSE
    )
  }
  
  vertical_positions <- trimws(
    as.character(vertical_positions)
  )
  
  vertical_positions <- unique(
    vertical_positions[
      nzchar(vertical_positions)
    ]
  )
  
  if (length(vertical_positions) == 0L) {
    stop(
      "At least one `vertical_positions` value is required.",
      call. = FALSE
    )
  }
  
  output_root <- normalizePath(
    output_root,
    mustWork = FALSE
  )
  
  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  combined_output_file <- file.path(
    output_root,
    "all_permafrost_all_depths_summary.csv"
  )
  
  status_output_file <- file.path(
    output_root,
    "all_permafrost_all_depths_status.csv"
  )
  
  message(
    "Permafrost lookup rows: ",
    nrow(lookup_dt)
  )
  message(
    "Vertical positions: ",
    paste(
      vertical_positions,
      collapse = ", "
    )
  )
  message(
    "Output root: ",
    output_root
  )
  
  results <- run_all_permafrost_indices_depths_parallel(
    lookup = lookup_dt,
    workers = workers,
    vertical_positions = vertical_positions,
    output_root = output_root,
    combined_output_file = combined_output_file,
    status_output_file = status_output_file,
    obs_start_year = as.integer(obs_start_year),
    obs_end_year = as.integer(obs_end_year),
    clim_root = clim_root,
    clim_basename = clim_basename,
    neon_env_dir = neon_env_dir,
    ...
  )
  
  saveRDS(
    results,
    file.path(
      output_root,
      "permafrost_tsoil_results.rds"
    )
  )
  
  if (file.exists(status_output_file)) {
    status <- data.table::fread(
      status_output_file
    )
    
    message(
      "Successful jobs: ",
      sum(
        status$status == "success",
        na.rm = TRUE
      ),
      " / ",
      nrow(status)
    )
  }
  
  message(
    "Combined summary: ",
    combined_output_file
  )
  message(
    "Combined status: ",
    status_output_file
  )
  
  invisible(
    results
  )
}


# =============================================================================
# Command-line entry point
# =============================================================================

run_permafrost_calibration_cli <- function(
    args = commandArgs(
      trailingOnly = TRUE
    )
) {
  if (length(args) < 1L || args[1L] %in% c("-h", "--help")) {
    stop(
      paste0(
        "Usage:\n",
        "  Rscript run_permafrost_calibration.R ",
        "<permafrost_lookup.csv|rds|RData> ",
        "[output_root] [workers] [depths]\n\n",
        "Example:\n",
        "  Rscript run_permafrost_calibration.R ",
        "/projectnb/dietzelab/guYANG/soilparam/permafrost_lookup.csv ",
        "/projectnb/dietzelab/guYANG/soilparam/permafrost_tsoil ",
        "16 501,502,503,504"
      ),
      call. = FALSE
    )
  }
  
  lookup_file <- args[1L]
  
  output_root <- if (length(args) >= 2L) {
    args[2L]
  } else {
    "/projectnb/dietzelab/guYANG/soilparam/permafrost_tsoil"
  }
  
  workers <- if (length(args) >= 3L) {
    suppressWarnings(
      as.integer(args[3L])
    )
  } else {
    16L
  }
  
  vertical_positions <- if (length(args) >= 4L) {
    strsplit(
      args[4L],
      split = ",",
      fixed = TRUE
    )[[1L]]
  } else {
    c(
      "501",
      "502",
      "503",
      "504"
    )
  }
  
  permafrost_lookup <- read_permafrost_lookup(
    lookup_file
  )
  
  run_permafrost_calibration(
    permafrost_lookup = permafrost_lookup,
    output_root = output_root,
    workers = workers,
    vertical_positions = vertical_positions
  )
}


if (sys.nframe() == 0L) {
  run_permafrost_calibration_cli()
}