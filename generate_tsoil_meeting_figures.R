#!/usr/bin/env Rscript

# =============================================================================
# Generate presentation-ready figures for the NEON SoilT / tau workflow
# =============================================================================
#
# All plot text is English. The default target is NEON vertical position 502
# (nominal depth approximately 6 cm).
#
# The script compares models on identical held-out observations:
#
# * Non-permafrost:
#     Original SIPNET SoilT and optimized-tau SoilT are read from each site's
#     LOYO validation file.
# * Permafrost:
#     Optimized predictions are read from the tau-only LOYO prediction table.
#     Original SIPNET SoilT is matched from the unmodified member-1 .clim file
#     at the same timestamps.
#
# Soil texture is read from the site-specific NetCDF file. For a requested
# point depth, the selected soil layer is the first layer whose bottom depth is
# greater than or equal to that point depth.
#
# Example:
#
# Rscript generate_tsoil_meeting_figures.R \
#   --depth=502 \
#   --output-dir=/projectnb/dietzelab/guYANG/Soil_Param/tsoil_meeting_figures_VER502
#
# @author Yang Gu


# =============================================================================
# Command-line configuration
# =============================================================================

parse_named_arguments <- function(arguments) {
  parsed <- list()
  
  for (argument in arguments) {
    if (!grepl("^--[^=]+=", argument)) {
      stop(
        "Every argument must use --name=value syntax. Invalid argument: ",
        argument,
        call. = FALSE
      )
    }
    
    key <- sub("^--([^=]+)=.*$", "\\1", argument)
    value <- sub("^--[^=]+=", "", argument)
    key <- gsub("-", "_", key, fixed = TRUE)
    parsed[[key]] <- value
  }
  
  parsed
}


arguments <- parse_named_arguments(
  commandArgs(trailingOnly = TRUE)
)


configuration <- list(
  depth = "502",
  lookup = "/projectnb/dietzelab/guYANG/soilparam/lookup.csv",
  nonpermafrost_summary = paste0(
    "/projectnb/dietzelab/guYANG/soilparam/tau_by_index/",
    "all_indices_all_depths_summary.csv"
  ),
  nonpermafrost_root =
    "/projectnb/dietzelab/guYANG/soilparam/tau_by_index",
  permafrost_summary = paste0(
    "/projectnb/dietzelab/guYANG/soilparam/permafrost_tau_only/",
    "permafrost_tsoil_tau_only_MLE_summary.csv"
  ),
  permafrost_predictions = paste0(
    "/projectnb/dietzelab/guYANG/soilparam/permafrost_tau_only/",
    "permafrost_tsoil_tau_only_LOYO_predictions.csv"
  ),
  clim_root =
    "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/ERA5_2012_2024",
  clim_basename = "ERA5.1.2012-01-01.2024-12-31.clim",
  soil_texture_root = paste0(
    "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/soil_nc/",
    "soil_texture_output/soil_texture_ensemble"
  ),
  output_dir =
    "/projectnb/dietzelab/guYANG/Soil_Param/tsoil_meeting_figures_VER502",
  match_tolerance_seconds = "60"
)


for (argument_name in names(arguments)) {
  if (!argument_name %in% names(configuration)) {
    stop(
      "Unknown argument --",
      gsub("_", "-", argument_name, fixed = TRUE),
      call. = FALSE
    )
  }
  
  configuration[[argument_name]] <- arguments[[argument_name]]
}


target_vertical_position <- as.character(configuration$depth)
target_depth_cm <- c(
  "501" = 2,
  "502" = 6,
  "503" = 16,
  "504" = 26
)[target_vertical_position]


if (length(target_depth_cm) != 1L || is.na(target_depth_cm)) {
  stop(
    "No nominal depth is defined for vertical position ",
    target_vertical_position,
    ". Supported values are 501, 502, 503, and 504.",
    call. = FALSE
  )
}


target_depth_cm <- as.numeric(target_depth_cm)
match_tolerance_seconds <- suppressWarnings(
  as.numeric(configuration$match_tolerance_seconds)
)


if (!is.finite(match_tolerance_seconds) || match_tolerance_seconds < 0) {
  stop("`match_tolerance_seconds` must be non-negative.", call. = FALSE)
}


# =============================================================================
# Dependencies and input validation
# =============================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "maps",
  "ncdf4"
)


missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1L),
    quietly = TRUE
  )
]


if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before rerunning this script.",
    call. = FALSE
  )
}


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


required_files <- c(
  lookup = configuration$lookup,
  nonpermafrost_summary = configuration$nonpermafrost_summary,
  permafrost_summary = configuration$permafrost_summary,
  permafrost_predictions = configuration$permafrost_predictions
)


missing_files <- required_files[!file.exists(required_files)]


if (length(missing_files) > 0L) {
  stop(
    "Required input files are missing:\n",
    paste(names(missing_files), missing_files, sep = ": ", collapse = "\n"),
    call. = FALSE
  )
}


dir.create(
  configuration$output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


message("Target vertical position: ", target_vertical_position)
message("Nominal depth: ", target_depth_cm, " cm")
message("Output directory: ", configuration$output_dir)


# =============================================================================
# General helpers
# =============================================================================

first_existing_column <- function(table, candidates, required = TRUE) {
  selected <- intersect(candidates, names(table))
  
  if (length(selected) == 0L) {
    if (isTRUE(required)) {
      stop(
        "None of the expected columns were found: ",
        paste(candidates, collapse = ", "),
        call. = FALSE
      )
    }
    
    return(NA_character_)
  }
  
  selected[1L]
}


filter_vertical_position <- function(
    table,
    vertical_position
) {
  table <- data.table::as.data.table(
    data.table::copy(
      table
    )
  )
  
  # Store the function argument under a name that cannot collide
  # with the data.table column named vertical_position.
  target_vertical_position_value <- as.character(
    vertical_position
  )[1L]
  
  if (
    is.na(target_vertical_position_value) ||
    !nzchar(target_vertical_position_value)
  ) {
    stop(
      "`vertical_position` must be one non-empty value.",
      call. = FALSE
    )
  }
  
  depth_column <- first_existing_column(
    table,
    c(
      "vertical_position",
      "verticalPosition"
    )
  )
  
  filtered_table <- table[
    as.character(
      get(depth_column)
    ) ==
      target_vertical_position_value
  ]
  
  if (nrow(filtered_table) == 0L) {
    stop(
      "No rows were found for vertical position ",
      target_vertical_position_value,
      ". Available positions: ",
      paste(
        sort(
          unique(
            as.character(
              table[[depth_column]]
            )
          )
        ),
        collapse = ", "
      ),
      call. = FALSE
    )
  }
  
  filtered_table
}

predictive_metrics <- function(observed, predicted) {
  valid <- is.finite(observed) & is.finite(predicted)
  observed <- as.numeric(observed[valid])
  predicted <- as.numeric(predicted[valid])
  n <- length(observed)
  
  if (n < 2L) {
    return(
      list(
        n = n,
        r2 = NA_real_,
        rmse_C = NA_real_,
        mae_C = NA_real_,
        bias_C = NA_real_,
        correlation = NA_real_
      )
    )
  }
  
  residual <- predicted - observed
  denominator <- sum((observed - mean(observed))^2)
  
  r2 <- if (is.finite(denominator) && denominator > 0) {
    1 - sum(residual^2) / denominator
  } else {
    NA_real_
  }
  
  correlation <- if (
    stats::sd(observed) > 0 &&
    stats::sd(predicted) > 0
  ) {
    stats::cor(observed, predicted)
  } else {
    NA_real_
  }
  
  list(
    n = n,
    r2 = r2,
    rmse_C = sqrt(mean(residual^2)),
    mae_C = mean(abs(residual)),
    bias_C = mean(residual),
    correlation = correlation
  )
}


compare_two_predictions <- function(
    table,
    group_columns,
    observation_column,
    original_column,
    improved_column
) {
  table <- data.table::as.data.table(data.table::copy(table))
  
  table <- table[
    is.finite(get(observation_column)) &
      is.finite(get(original_column)) &
      is.finite(get(improved_column))
  ]
  
  table[
    ,
    {
      original_metrics <- predictive_metrics(
        get(observation_column),
        get(original_column)
      )
      
      improved_metrics <- predictive_metrics(
        get(observation_column),
        get(improved_column)
      )
      
      list(
        n_validation = original_metrics$n,
        original_r2 = original_metrics$r2,
        improved_r2 = improved_metrics$r2,
        delta_r2 = improved_metrics$r2 - original_metrics$r2,
        original_rmse_C = original_metrics$rmse_C,
        improved_rmse_C = improved_metrics$rmse_C,
        delta_rmse_C = improved_metrics$rmse_C - original_metrics$rmse_C,
        original_mae_C = original_metrics$mae_C,
        improved_mae_C = improved_metrics$mae_C,
        original_bias_C = original_metrics$bias_C,
        improved_bias_C = improved_metrics$bias_C,
        original_correlation = original_metrics$correlation,
        improved_correlation = improved_metrics$correlation
      )
    },
    by = group_columns
  ]
}


parse_utc_time <- function(values) {
  if (inherits(values, "POSIXct")) {
    return(as.POSIXct(values, tz = "UTC"))
  }
  
  if (is.numeric(values)) {
    return(as.POSIXct(values, origin = "1970-01-01", tz = "UTC"))
  }
  
  as.POSIXct(as.character(values), tz = "UTC")
}


nearest_match_values <- function(
    target_time,
    reference_time,
    reference_value,
    tolerance_seconds
) {
  target_numeric <- as.numeric(parse_utc_time(target_time))
  reference_numeric <- as.numeric(parse_utc_time(reference_time))
  
  valid_target <- is.finite(target_numeric)
  result <- rep(NA_real_, length(target_numeric))
  
  if (!any(valid_target)) {
    return(result)
  }
  
  valid_reference <- is.finite(reference_numeric) & is.finite(reference_value)
  reference_numeric <- reference_numeric[valid_reference]
  reference_value <- as.numeric(reference_value[valid_reference])
  
  if (length(reference_numeric) == 0L) {
    return(result)
  }
  
  order_index <- order(reference_numeric)
  reference_numeric <- reference_numeric[order_index]
  reference_value <- reference_value[order_index]
  
  target_valid_numeric <- target_numeric[valid_target]
  interval <- findInterval(target_valid_numeric, reference_numeric)
  left <- pmax(interval, 1L)
  right <- pmin(interval + 1L, length(reference_numeric))
  
  left_distance <- abs(target_valid_numeric - reference_numeric[left])
  right_distance <- abs(target_valid_numeric - reference_numeric[right])
  use_right <- right_distance < left_distance
  selected <- left
  selected[use_right] <- right[use_right]
  
  distance <- abs(target_valid_numeric - reference_numeric[selected])
  matched <- reference_value[selected]
  matched[!is.finite(distance) | distance > tolerance_seconds] <- NA_real_
  result[valid_target] <- matched
  result
}


read_sipnet_clim_minimal <- function(path) {
  if (!file.exists(path)) {
    stop("SIPNET climate file not found: ", path, call. = FALSE)
  }
  
  clim <- data.table::fread(
    path,
    header = FALSE,
    showProgress = FALSE
  )
  
  names_12 <- c(
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
  
  if (ncol(clim) == 12L) {
    data.table::setnames(clim, names_12)
  } else if (ncol(clim) == 14L) {
    data.table::setnames(
      clim,
      c("grid_index", names_12, "soil_wetness")
    )
  } else {
    stop(
      "Expected 12 or 14 columns in ",
      path,
      "; found ",
      ncol(clim),
      ".",
      call. = FALSE
    )
  }
  
  clim[
    ,
    time := as.POSIXct(
      sprintf("%04d-01-01 00:00:00", as.integer(year)),
      tz = "UTC"
    ) +
      (as.numeric(doy) - 1) * 86400 +
      as.numeric(hour) * 3600
  ]
  
  clim[
    ,
    .(
      time,
      Tsoil_original_C = as.numeric(soil_temp_C)
    )
  ]
}


meeting_theme <- function(base_size = 15) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 3,
        margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey30",
        margin = ggplot2::margin(b = 12)
      ),
      plot.caption = ggplot2::element_text(
        color = "grey40",
        hjust = 0
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold")
    )
}


save_figure <- function(plot, stem, width, height) {
  png_file <- file.path(
    configuration$output_dir,
    paste0(stem, ".png")
  )
  
  pdf_file <- file.path(
    configuration$output_dir,
    paste0(stem, ".pdf")
  )
  
  ggplot2::ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )
  
  tryCatch(
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    ),
    error = function(error) {
      warning(
        "Cairo PDF output failed; using the standard PDF device: ",
        conditionMessage(error)
      )
      
      ggplot2::ggsave(
        filename = pdf_file,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = grDevices::pdf,
        bg = "white"
      )
    }
  )
  
  data.table::data.table(
    figure = stem,
    png_file = png_file,
    pdf_file = pdf_file
  )
}


add_site_labels <- function(plot, data, label_column = "NEON_code") {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    plot +
      ggrepel::geom_text_repel(
        data = data,
        ggplot2::aes(label = get(label_column)),
        size = 3.4,
        min.segment.length = 0,
        max.overlaps = Inf,
        box.padding = 0.35,
        point.padding = 0.2,
        seed = 341
      )
  } else {
    plot +
      ggplot2::geom_text(
        data = data,
        ggplot2::aes(label = get(label_column)),
        size = 3,
        nudge_y = 0.6,
        check_overlap = TRUE
      )
  }
}


# =============================================================================
# Lookup and parameter summaries
# =============================================================================

lookup <- data.table::fread(configuration$lookup)

latitude_column <- first_existing_column(
  lookup,
  c("site_lat", "latitude", "lat", "LAT")
)

longitude_column <- first_existing_column(
  lookup,
  c("site_lon", "longitude", "lon", "LONG", "LON")
)

site_column <- first_existing_column(
  lookup,
  c("NEON_code", "site", "siteID")
)


lookup_standard <- lookup[
  ,
  .(
    index = as.integer(index),
    NEON_code = as.character(get(site_column)),
    latitude = as.numeric(get(latitude_column)),
    longitude = as.numeric(get(longitude_column))
  )
]


lookup_standard <- unique(lookup_standard, by = "index")


non_summary <- filter_vertical_position(
  data.table::fread(configuration$nonpermafrost_summary),
  target_vertical_position
)


permafrost_summary <- filter_vertical_position(
  data.table::fread(configuration$permafrost_summary),
  target_vertical_position
)


non_summary[
  ,
  `:=`(
    index = as.integer(index),
    tau_days = as.numeric(tau_days)
  )
]


permafrost_summary[
  ,
  `:=`(
    index = as.integer(index),
    tau_days = as.numeric(tau_days)
  )
]


if (anyDuplicated(non_summary$index) > 0L) {
  stop(
    "The non-permafrost summary has duplicate index rows at VER",
    target_vertical_position,
    ".",
    call. = FALSE
  )
}


if (anyDuplicated(permafrost_summary$index) > 0L) {
  stop(
    "The permafrost summary has duplicate index rows at VER",
    target_vertical_position,
    ".",
    call. = FALSE
  )
}


# =============================================================================
# Non-permafrost held-out comparison
# =============================================================================

non_validation_files <- list.files(
  configuration$nonpermafrost_root,
  pattern = "_soil_temperature_validation\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)


read_nonpermafrost_validation <- function(path) {
  directory_name <- basename(dirname(path))
  
  if (!grepl("_index[0-9]+_VER[^/]+$", directory_name)) {
    return(NULL)
  }
  
  index_value <- suppressWarnings(
    as.integer(sub(".*_index([0-9]+)_VER.*$", "\\1", directory_name))
  )
  
  depth_value <- sub(".*_VER([^/]+)$", "\\1", directory_name)
  site_value <- sub("_index[0-9]+_VER.*$", "", directory_name)
  
  if (
    !is.finite(index_value) ||
    depth_value != target_vertical_position
  ) {
    return(NULL)
  }
  
  validation <- tryCatch(
    data.table::fread(path),
    error = function(error) {
      warning("Could not read ", path, ": ", conditionMessage(error))
      NULL
    }
  )
  
  if (is.null(validation) || nrow(validation) == 0L) {
    return(NULL)
  }
  
  required_columns <- c(
    "Tsoil_obs_C",
    "Tsoil_current_C",
    "Tsoil_tau_opt_C"
  )
  
  if (!all(required_columns %in% names(validation))) {
    warning(
      "Skipping ",
      path,
      "; required LOYO prediction columns are missing."
    )
    return(NULL)
  }
  
  validation[
    ,
    `:=`(
      index = index_value,
      NEON_code = site_value
    )
  ]
  
  validation
}


non_validation <- data.table::rbindlist(
  lapply(non_validation_files, read_nonpermafrost_validation),
  use.names = TRUE,
  fill = TRUE
)


if (nrow(non_validation) == 0L) {
  stop(
    "No usable non-permafrost LOYO validation files were found under ",
    configuration$nonpermafrost_root,
    ".",
    call. = FALSE
  )
}


non_metrics <- compare_two_predictions(
  table = non_validation,
  group_columns = c("index", "NEON_code"),
  observation_column = "Tsoil_obs_C",
  original_column = "Tsoil_current_C",
  improved_column = "Tsoil_tau_opt_C"
)


non_sites <- merge(
  non_summary[
    ,
    .(
      index,
      tau_days,
      tau_ci_low = as.numeric(tau_ci_low_profile_approx),
      tau_ci_high = as.numeric(tau_ci_high_profile_approx)
    )
  ],
  non_metrics,
  by = "index",
  all.x = TRUE
)


non_sites <- merge(
  non_sites,
  lookup_standard,
  by = "index",
  all.x = TRUE,
  suffixes = c("", "_lookup")
)


if ("NEON_code_lookup" %in% names(non_sites)) {
  non_sites[
    is.na(NEON_code) | !nzchar(NEON_code),
    NEON_code := NEON_code_lookup
  ]
  non_sites[, NEON_code_lookup := NULL]
}


non_sites[
  ,
  `:=`(
    model_group = "Non-permafrost",
    validation_resolution = "Daily"
  )
]


# =============================================================================
# Permafrost held-out comparison on identical 3-hour observations
# =============================================================================

permafrost_predictions <- filter_vertical_position(
  data.table::fread(configuration$permafrost_predictions),
  target_vertical_position
)


required_permafrost_columns <- c(
  "index",
  "time",
  "Tsoil_obs_C",
  "Tsoil_pred_C"
)


if (!all(required_permafrost_columns %in% names(permafrost_predictions))) {
  stop(
    "Permafrost LOYO predictions are missing: ",
    paste(
      setdiff(required_permafrost_columns, names(permafrost_predictions)),
      collapse = ", "
    ),
    call. = FALSE
  )
}


permafrost_predictions[
  ,
  `:=`(
    index = as.integer(index),
    time = parse_utc_time(time),
    Tsoil_obs_C = as.numeric(Tsoil_obs_C),
    Tsoil_pred_C = as.numeric(Tsoil_pred_C)
  )
]


permafrost_predictions[, Tsoil_original_C := NA_real_]


for (site_index in sort(unique(permafrost_predictions$index))) {
  clim_file <- file.path(
    configuration$clim_root,
    sprintf("ERA5_%d_1", site_index),
    configuration$clim_basename
  )
  
  original_clim <- read_sipnet_clim_minimal(clim_file)
  rows <- which(permafrost_predictions$index == site_index)
  
  permafrost_predictions[
    rows,
    Tsoil_original_C := nearest_match_values(
      target_time = time,
      reference_time = original_clim$time,
      reference_value = original_clim$Tsoil_original_C,
      tolerance_seconds = match_tolerance_seconds
    )
  ]
}


permafrost_predictions <- merge(
  permafrost_predictions,
  lookup_standard[, .(index, NEON_code)],
  by = "index",
  all.x = TRUE
)


permafrost_metrics <- compare_two_predictions(
  table = permafrost_predictions,
  group_columns = c("index", "NEON_code"),
  observation_column = "Tsoil_obs_C",
  original_column = "Tsoil_original_C",
  improved_column = "Tsoil_pred_C"
)


permafrost_sites <- merge(
  permafrost_summary[
    ,
    .(
      index,
      tau_days,
      fixed_a_C = if ("fixed_a_C" %in% names(permafrost_summary)) {
        as.numeric(fixed_a_C)
      } else {
        NA_real_
      },
      fixed_n_warm = if ("fixed_n_warm" %in% names(permafrost_summary)) {
        as.numeric(fixed_n_warm)
      } else {
        NA_real_
      },
      fixed_n_cold = if ("fixed_n_cold" %in% names(permafrost_summary)) {
        as.numeric(fixed_n_cold)
      } else {
        NA_real_
      }
    )
  ],
  permafrost_metrics,
  by = "index",
  all.x = TRUE
)


permafrost_sites <- merge(
  permafrost_sites,
  lookup_standard,
  by = "index",
  all.x = TRUE,
  suffixes = c("", "_lookup")
)


if ("NEON_code_lookup" %in% names(permafrost_sites)) {
  permafrost_sites[
    is.na(NEON_code) | !nzchar(NEON_code),
    NEON_code := NEON_code_lookup
  ]
  permafrost_sites[, NEON_code_lookup := NULL]
}


permafrost_sites[
  ,
  `:=`(
    model_group = "Permafrost",
    validation_resolution = "3-hour"
  )
]


# =============================================================================
# Combined site-level metrics
# =============================================================================

site_metrics <- data.table::rbindlist(
  list(non_sites, permafrost_sites),
  use.names = TRUE,
  fill = TRUE
)


data.table::setcolorder(
  site_metrics,
  c(
    "model_group",
    "index",
    "NEON_code",
    "latitude",
    "longitude",
    "tau_days",
    "n_validation",
    "original_r2",
    "improved_r2",
    "delta_r2",
    setdiff(
      names(site_metrics),
      c(
        "model_group",
        "index",
        "NEON_code",
        "latitude",
        "longitude",
        "tau_days",
        "n_validation",
        "original_r2",
        "improved_r2",
        "delta_r2"
      )
    )
  )
)


site_metrics_file <- file.path(
  configuration$output_dir,
  paste0("site_tsoil_metrics_VER", target_vertical_position, ".csv")
)


data.table::fwrite(site_metrics, site_metrics_file)


# =============================================================================
# Soil texture at the target depth and non-permafrost regression
# =============================================================================

find_soil_texture_file <- function(site_index) {
  expected <- file.path(
    configuration$soil_texture_root,
    as.character(site_index),
    sprintf("Soil_params_0-%d_1.nc", site_index)
  )
  
  if (file.exists(expected)) {
    return(expected)
  }
  
  directory <- file.path(
    configuration$soil_texture_root,
    as.character(site_index)
  )
  
  candidates <- if (dir.exists(directory)) {
    list.files(
      directory,
      pattern = sprintf("-%d_1\\.nc$", site_index),
      full.names = TRUE
    )
  } else {
    character()
  }
  
  if (length(candidates) == 0L) {
    return(NA_character_)
  }
  
  candidates[1L]
}


read_soil_texture <- function(site_index, point_depth_cm) {
  path <- find_soil_texture_file(site_index)
  
  if (is.na(path) || !file.exists(path)) {
    warning("No soil texture NetCDF found for index ", site_index, ".")
    return(
      data.table::data.table(
        index = as.integer(site_index),
        soil_texture_file = NA_character_,
        requested_depth_cm = point_depth_cm,
        selected_layer = NA_integer_,
        layer_bottom_depth_m = NA_real_,
        sand_pct = NA_real_,
        clay_pct = NA_real_,
        silt_pct = NA_real_
      )
    )
  }
  
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  depth_values <- nc$dim$depth$vals
  point_depth_m <- point_depth_cm / 100
  
  selected_layer <- which(depth_values >= point_depth_m)[1L]
  
  if (is.na(selected_layer)) {
    selected_layer <- which.min(abs(depth_values - point_depth_m))
  }
  
  read_fraction <- function(variable_name) {
    values <- ncdf4::ncvar_get(nc, variable_name)
    value <- as.numeric(values[selected_layer])
    
    if (is.finite(value) && abs(value) <= 1.5) {
      value <- value * 100
    }
    
    value
  }
  
  data.table::data.table(
    index = as.integer(site_index),
    soil_texture_file = path,
    requested_depth_cm = point_depth_cm,
    selected_layer = as.integer(selected_layer),
    layer_bottom_depth_m = as.numeric(depth_values[selected_layer]),
    sand_pct = read_fraction("fraction_of_sand_in_soil"),
    clay_pct = read_fraction("fraction_of_clay_in_soil"),
    silt_pct = read_fraction("fraction_of_silt_in_soil")
  )
}


soil_texture <- data.table::rbindlist(
  lapply(
    non_sites$index,
    read_soil_texture,
    point_depth_cm = target_depth_cm
  ),
  use.names = TRUE,
  fill = TRUE
)


non_texture <- merge(
  non_sites,
  soil_texture,
  by = "index",
  all.x = TRUE
)


non_texture_file <- file.path(
  configuration$output_dir,
  paste0(
    "nonpermafrost_tau_soil_texture_VER",
    target_vertical_position,
    ".csv"
  )
)


data.table::fwrite(non_texture, non_texture_file)


regression_data <- non_texture[
  is.finite(tau_days) &
    tau_days > 0 &
    is.finite(sand_pct) &
    is.finite(clay_pct) &
    is.finite(silt_pct)
]


if (nrow(regression_data) < 5L) {
  stop(
    "Fewer than five non-permafrost sites have complete tau and soil texture data.",
    call. = FALSE
  )
}


# Sand + clay + silt sum to approximately 100%. Silt is the reference
# component and is omitted from the multiple regression to avoid perfect
# compositional collinearity.
raw_texture_model <- stats::lm(
  tau_days ~ sand_pct + clay_pct,
  data = regression_data
)


log_texture_model <- stats::lm(
  log(tau_days) ~ sand_pct + clay_pct,
  data = regression_data
)


regression_data[
  ,
  `:=`(
    tau_fitted_raw_days = as.numeric(stats::predict(raw_texture_model)),
    tau_fitted_log_days = exp(as.numeric(stats::predict(log_texture_model)))
  )
]


raw_coefficients <- stats::coef(raw_texture_model)
log_coefficients <- stats::coef(log_texture_model)
raw_model_summary <- summary(raw_texture_model)
log_model_summary <- summary(log_texture_model)


signed_term <- function(value, digits = 3) {
  sprintf(
    if (value >= 0) "+ %.*f" else "- %.*f",
    digits,
    abs(value)
  )
}


raw_equation <- paste0(
  "Tau = ",
  sprintf("%.2f", raw_coefficients[1L]),
  " ",
  signed_term(raw_coefficients[2L]),
  " × Sand (%) ",
  signed_term(raw_coefficients[3L]),
  " × Clay (%)"
)


log_equation <- paste0(
  "log(Tau) = ",
  sprintf("%.3f", log_coefficients[1L]),
  " ",
  signed_term(log_coefficients[2L], digits = 4),
  " × Sand (%) ",
  signed_term(log_coefficients[3L], digits = 4),
  " × Clay (%)"
)


coefficient_table <- data.table::rbindlist(
  list(
    data.table::data.table(
      model = "Raw tau",
      term = names(raw_coefficients),
      estimate = as.numeric(raw_coefficients)
    ),
    data.table::data.table(
      model = "Log tau",
      term = names(log_coefficients),
      estimate = as.numeric(log_coefficients)
    )
  )
)


data.table::fwrite(
  coefficient_table,
  file.path(
    configuration$output_dir,
    paste0(
      "nonpermafrost_tau_texture_coefficients_VER",
      target_vertical_position,
      ".csv"
    )
  )
)


regression_report <- c(
  "NON-PERMAFROST TAU ~ SOIL TEXTURE REGRESSION",
  paste0("NEON vertical position: ", target_vertical_position),
  paste0("Nominal sensor depth: ", target_depth_cm, " cm"),
  paste0("Complete sites: ", nrow(regression_data)),
  "",
  "Primary raw-scale multiple regression",
  raw_equation,
  sprintf("R-squared = %.4f", raw_model_summary$r.squared),
  sprintf("Adjusted R-squared = %.4f", raw_model_summary$adj.r.squared),
  "",
  "Log-scale sensitivity regression",
  log_equation,
  sprintf("R-squared = %.4f", log_model_summary$r.squared),
  sprintf("Adjusted R-squared = %.4f", log_model_summary$adj.r.squared),
  "",
  paste0(
    "Interpretation note: sand + clay + silt sum to approximately 100%; ",
    "silt is omitted and acts as the reference component to avoid perfect ",
    "collinearity. Coefficients are per one percentage-point change in texture."
  ),
  paste0(
    "Layer-selection rule: first NetCDF layer whose bottom depth is at or ",
    "below the requested point depth."
  )
)


writeLines(
  regression_report,
  file.path(
    configuration$output_dir,
    paste0(
      "nonpermafrost_tau_soil_texture_regression_VER",
      target_vertical_position,
      ".txt"
    )
  )
)


# =============================================================================
# Figures
# =============================================================================

figure_manifest <- list()


world_map <- ggplot2::map_data("world")
north_america <- world_map[
  world_map$region %in% c("USA", "Canada", "Mexico"),
]


plot_delta_r2_map <- function(data, group_title) {
  plot_data <- data[
    is.finite(longitude) &
      is.finite(latitude) &
      is.finite(delta_r2)
  ]
  
  if (nrow(plot_data) == 0L) {
    stop("No mappable delta R-squared values for ", group_title, ".")
  }
  
  symmetric_limit <- max(abs(plot_data$delta_r2), na.rm = TRUE)
  symmetric_limit <- max(symmetric_limit, 0.05)
  resolution_label <- unique(plot_data$validation_resolution)
  resolution_label <- resolution_label[!is.na(resolution_label)][1L]
  
  plot <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = north_america,
      ggplot2::aes(x = long, y = lat, group = group),
      fill = "grey94",
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(
        x = longitude,
        y = latitude,
        color = delta_r2
      ),
      size = 4.5,
      alpha = 0.95
    ) +
    ggplot2::scale_color_gradient2(
      low = "#B2182B",
      mid = "#F7F7F7",
      high = "#2166AC",
      midpoint = 0,
      limits = c(-symmetric_limit, symmetric_limit),
      name = expression(Delta * R^2)
    ) +
    ggplot2::coord_quickmap(
      xlim = c(-170, -52),
      ylim = c(20, 75),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = paste0(group_title, ": change in held-out R²"),
      subtitle = paste0(
        "Optimized site-specific tau minus original SIPNET SoilT | NEON VER",
        target_vertical_position,
        " (~",
        target_depth_cm,
        " cm) | ",
        resolution_label,
        " held-out observations"
      ),
      x = "Longitude",
      y = "Latitude",
      caption = paste0(
        "Positive values indicate improved out-of-sample performance. ",
        "R² is calculated on identical held-out observations."
      )
    ) +
    meeting_theme(15) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.key.width = grid::unit(2.2, "cm")
    )
  
  add_site_labels(plot, plot_data)
}


non_map <- plot_delta_r2_map(non_sites, "Non-permafrost sites")
figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  non_map,
  paste0("01_nonpermafrost_delta_R2_map_VER", target_vertical_position),
  width = 13,
  height = 8.5
)


permafrost_map <- plot_delta_r2_map(permafrost_sites, "Permafrost sites")
figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  permafrost_map,
  paste0("02_permafrost_delta_R2_map_VER", target_vertical_position),
  width = 13,
  height = 8.5
)


plot_tau_by_site <- function(data, group_title, point_color) {
  plot_data <- data[
    is.finite(tau_days) & tau_days > 0
  ]
  
  plot_data[
    ,
    site_label := sprintf("%s (%d)", NEON_code, index)
  ]
  
  plot_data[
    ,
    site_label := factor(
      site_label,
      levels = site_label[order(tau_days)]
    )
  ]
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = tau_days, y = site_label)
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = min(tau_days) * 0.8, xend = tau_days, yend = site_label),
      color = "grey75",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      color = point_color,
      size = 3.8
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f d", tau_days)),
      hjust = -0.18,
      size = 3.4,
      color = "grey20"
    ) +
    ggplot2::scale_x_log10(
      expand = ggplot2::expansion(mult = c(0.04, 0.22))
    ) +
    ggplot2::labs(
      title = paste0(group_title, ": site-specific thermal memory"),
      subtitle = paste0(
        "Final MLE tau at NEON VER",
        target_vertical_position,
        " (~",
        target_depth_cm,
        " cm)"
      ),
      x = "Tau (days; logarithmic scale)",
      y = "NEON site (model index)",
      caption = "Larger tau indicates a slower and smoother soil-temperature response to air temperature."
    ) +
    meeting_theme(14)
}


non_tau_plot <- plot_tau_by_site(
  non_sites,
  "Non-permafrost sites",
  "#0072B2"
)


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  non_tau_plot,
  paste0("03_nonpermafrost_tau_by_site_VER", target_vertical_position),
  width = 11,
  height = max(8, 4 + 0.22 * nrow(non_sites))
)


permafrost_tau_plot <- plot_tau_by_site(
  permafrost_sites,
  "Permafrost sites",
  "#D55E00"
)


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  permafrost_tau_plot,
  paste0("04_permafrost_tau_by_site_VER", target_vertical_position),
  width = 10,
  height = max(5, 3 + 0.45 * nrow(permafrost_sites))
)


r2_plot_data <- site_metrics[
  is.finite(original_r2) & is.finite(improved_r2)
]


r2_range <- range(
  c(r2_plot_data$original_r2, r2_plot_data$improved_r2),
  na.rm = TRUE
)


r2_padding <- max(diff(r2_range) * 0.08, 0.05)


r2_comparison_plot <- ggplot2::ggplot(
  r2_plot_data,
  ggplot2::aes(
    x = original_r2,
    y = improved_r2,
    color = model_group
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey45"
  ) +
  ggplot2::geom_point(size = 3.8, alpha = 0.9) +
  ggplot2::facet_wrap(~model_group) +
  ggplot2::coord_equal(
    xlim = c(r2_range[1L] - r2_padding, r2_range[2L] + r2_padding),
    ylim = c(r2_range[1L] - r2_padding, r2_range[2L] + r2_padding)
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Non-permafrost" = "#0072B2",
      "Permafrost" = "#D55E00"
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Held-out SoilT performance: original versus optimized",
    subtitle = paste0(
      "Points above the 1:1 line improved after site-specific tau calibration | VER",
      target_vertical_position
    ),
    x = "Original SIPNET SoilT R²",
    y = "Optimized SoilT R²",
    caption = paste0(
      "Predictive R² = 1 − SSE/SST; values can be negative. ",
      "Non-permafrost validation is daily; permafrost validation is 3-hour."
    )
  ) +
  meeting_theme(15)


r2_label_data <- r2_plot_data[
  order(-abs(delta_r2)),
  head(.SD, min(.N, 7L)),
  by = model_group
]


if (requireNamespace("ggrepel", quietly = TRUE)) {
  r2_comparison_plot <- r2_comparison_plot +
    ggrepel::geom_text_repel(
      data = r2_label_data,
      ggplot2::aes(label = NEON_code),
      size = 3.3,
      max.overlaps = Inf,
      seed = 502,
      show.legend = FALSE
    )
}


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  r2_comparison_plot,
  paste0("05_original_vs_optimized_R2_VER", target_vertical_position),
  width = 12,
  height = 6.5
)


delta_tau_plot <- ggplot2::ggplot(
  site_metrics[
    is.finite(tau_days) &
      tau_days > 0 &
      is.finite(delta_r2)
  ],
  ggplot2::aes(
    x = tau_days,
    y = delta_r2,
    color = model_group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45"
  ) +
  ggplot2::geom_point(size = 4, alpha = 0.9) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_color_manual(
    values = c(
      "Non-permafrost" = "#0072B2",
      "Permafrost" = "#D55E00"
    ),
    name = "Model group"
  ) +
  ggplot2::labs(
    title = "Does the estimated thermal memory explain model improvement?",
    subtitle = paste0("Site-level results at NEON VER", target_vertical_position),
    x = "Tau (days; logarithmic scale)",
    y = expression(Delta * R^2 == R[new]^2 - R[original]^2),
    caption = "Each point is one NEON site; values above zero indicate improvement."
  ) +
  meeting_theme(15)


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  delta_tau_plot,
  paste0("06_delta_R2_vs_tau_VER", target_vertical_position),
  width = 11,
  height = 7
)


regression_annotation <- paste0(
  raw_equation,
  "\nR² = ",
  sprintf("%.3f", raw_model_summary$r.squared),
  "; adjusted R² = ",
  sprintf("%.3f", raw_model_summary$adj.r.squared),
  "; n = ",
  nrow(regression_data)
)


regression_plot <- ggplot2::ggplot(
  regression_data,
  ggplot2::aes(
    x = tau_fitted_raw_days,
    y = tau_days
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "grey45"
  ) +
  ggplot2::geom_point(
    color = "#0072B2",
    size = 4,
    alpha = 0.9
  ) +
  ggplot2::annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = regression_annotation,
    hjust = -0.03,
    vjust = 1.15,
    size = 4.2,
    label.size = 0.25,
    fill = "white"
  ) +
  ggplot2::labs(
    title = "Non-permafrost tau predicted from soil texture",
    subtitle = paste0(
      "Multiple linear regression at NEON VER",
      target_vertical_position,
      " (~",
      target_depth_cm,
      " cm); silt is the reference component"
    ),
    x = "Fitted tau from sand and clay (days)",
    y = "MLE tau (days)",
    caption = paste0(
      "Sand, clay, and silt are compositional; silt is omitted to avoid perfect collinearity. ",
      "The full equation is also written to a text file."
    )
  ) +
  meeting_theme(15)


if (requireNamespace("ggrepel", quietly = TRUE)) {
  regression_plot <- regression_plot +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = NEON_code),
      size = 3.1,
      max.overlaps = Inf,
      seed = 616
    )
}


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  regression_plot,
  paste0(
    "07_nonpermafrost_tau_soil_texture_regression_VER",
    target_vertical_position
  ),
  width = 11,
  height = 7.5
)


texture_long <- data.table::melt(
  regression_data[
    ,
    .(
      index,
      NEON_code,
      tau_days,
      sand_pct,
      clay_pct,
      silt_pct
    )
  ],
  id.vars = c("index", "NEON_code", "tau_days"),
  measure.vars = c("sand_pct", "clay_pct", "silt_pct"),
  variable.name = "texture_component",
  value.name = "texture_pct"
)


texture_long[
  ,
  texture_component := factor(
    texture_component,
    levels = c("sand_pct", "clay_pct", "silt_pct"),
    labels = c("Sand", "Clay", "Silt")
  )
]


component_statistics <- texture_long[
  ,
  {
    model <- stats::lm(tau_days ~ texture_pct)
    coefficients <- stats::coef(model)
    
    list(
      intercept = coefficients[1L],
      slope = coefficients[2L],
      r2 = summary(model)$r.squared,
      label = sprintf(
        "Tau = %.2f %+.3f × texture\nR² = %.3f",
        coefficients[1L],
        coefficients[2L],
        summary(model)$r.squared
      )
    )
  },
  by = texture_component
]


data.table::fwrite(
  component_statistics,
  file.path(
    configuration$output_dir,
    paste0(
      "nonpermafrost_tau_texture_component_regressions_VER",
      target_vertical_position,
      ".csv"
    )
  )
)


texture_component_plot <- ggplot2::ggplot(
  texture_long,
  ggplot2::aes(x = texture_pct, y = tau_days)
) +
  ggplot2::geom_point(
    color = "#0072B2",
    size = 3.2,
    alpha = 0.8
  ) +
  ggplot2::geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#D55E00",
    fill = "#F0E442",
    linewidth = 0.9,
    alpha = 0.25
  ) +
  ggplot2::geom_text(
    data = component_statistics,
    ggplot2::aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.2,
    size = 3.7
  ) +
  ggplot2::facet_wrap(~texture_component, scales = "free_x") +
  ggplot2::labs(
    title = "Non-permafrost tau versus individual soil texture components",
    subtitle = paste0(
      "Univariate diagnostic relationships at NEON VER",
      target_vertical_position,
      " (~",
      target_depth_cm,
      " cm)"
    ),
    x = "Soil texture component (%)",
    y = "MLE tau (days)",
    caption = "These univariate fits are descriptive; the primary multiple regression uses sand and clay together."
  ) +
  meeting_theme(14)


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  texture_component_plot,
  paste0(
    "08_nonpermafrost_tau_vs_texture_components_VER",
    target_vertical_position
  ),
  width = 14,
  height = 6.5
)


formula_plot <- ggplot2::ggplot() +
  ggplot2::xlim(0, 1) +
  ggplot2::ylim(0, 1) +
  ggplot2::annotate(
    "text",
    x = 0.03,
    y = 0.94,
    hjust = 0,
    label = "Soil temperature models and validation target",
    size = 7,
    fontface = "bold"
  ) +
  ggplot2::annotate(
    "label",
    x = 0.03,
    y = 0.78,
    hjust = 0,
    label = paste0(
      "Non-permafrost model\n",
      "Tsoil(t) = Tsoil(t−1) + [1 − exp(−Δt/τ)] × [Tair(t) − Tsoil(t−1)]"
    ),
    size = 5.1,
    lineheight = 1.25,
    fill = "#E6F2FA",
    color = "#004C6D",
    label.size = 0.3
  ) +
  ggplot2::annotate(
    "label",
    x = 0.03,
    y = 0.54,
    hjust = 0,
    label = paste0(
      "Permafrost model\n",
      "Teff(t) = a + nwarm × max[Tair(t), 0] + ncold × min[Tair(t), 0]\n",
      "Tsoil(t) = Tsoil(t−1) + [1 − exp(−Δt/τ)] × [Teff(t) − Tsoil(t−1)]"
    ),
    size = 4.8,
    lineheight = 1.2,
    fill = "#FCE8DE",
    color = "#8C2D04",
    label.size = 0.3
  ) +
  ggplot2::annotate(
    "label",
    x = 0.03,
    y = 0.29,
    hjust = 0,
    label = paste0(
      "Original benchmark\n",
      "Soil temperature stored in the unmodified PEcAn-generated SIPNET climate file"
    ),
    size = 4.7,
    lineheight = 1.2,
    fill = "grey95",
    color = "grey20",
    label.size = 0.3
  ) +
  ggplot2::annotate(
    "text",
    x = 0.03,
    y = 0.10,
    hjust = 0,
    label = "Validation metric:  ΔR² = R²(optimized, held out) − R²(original, held out)",
    size = 5.2,
    fontface = "bold",
    color = "#333333"
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = "white", color = NA)
  )


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  formula_plot,
  paste0("09_tsoil_model_equations_VER", target_vertical_position),
  width = 14,
  height = 8
)


# Representative held-out time series. For each model group, select the site
# whose delta R-squared is closest to the group median rather than selecting an
# unusually good or bad case.
choose_representative_index <- function(data) {
  candidates <- data[is.finite(delta_r2)]
  
  if (nrow(candidates) == 0L) {
    return(NA_integer_)
  }
  
  target <- stats::median(candidates$delta_r2)
  candidates$index[which.min(abs(candidates$delta_r2 - target))]
}


representative_non_index <- choose_representative_index(non_sites)
representative_permafrost_index <- choose_representative_index(
  permafrost_sites
)


non_representative <- non_validation[
  index == representative_non_index &
    is.finite(Tsoil_obs_C) &
    is.finite(Tsoil_current_C) &
    is.finite(Tsoil_tau_opt_C)
]


non_representative[, plot_date := as.Date(date)]


non_year <- non_representative[
  ,
  .N,
  by = held_out_year
][
  which.max(N),
  held_out_year
]


non_representative <- non_representative[
  held_out_year == non_year
]


non_delta <- non_sites[
  index == representative_non_index,
  delta_r2
][1L]


non_site_name <- non_sites[
  index == representative_non_index,
  NEON_code
][1L]


non_time_series <- non_representative[
  ,
  .(
    plot_date,
    `NEON observations` = Tsoil_obs_C,
    `Original SIPNET` = Tsoil_current_C,
    `Optimized tau` = Tsoil_tau_opt_C
  )
]


non_time_series[
  ,
  panel := sprintf(
    "Non-permafrost: %s | held-out %s | ΔR² = %.3f",
    non_site_name,
    non_year,
    non_delta
  )
]


permafrost_representative <- permafrost_predictions[
  index == representative_permafrost_index &
    is.finite(Tsoil_obs_C) &
    is.finite(Tsoil_original_C) &
    is.finite(Tsoil_pred_C)
]


permafrost_year <- permafrost_representative[
  ,
  .N,
  by = held_out_year
][
  which.max(N),
  held_out_year
]


permafrost_representative <- permafrost_representative[
  held_out_year == permafrost_year
]


# Aggregate the permafrost example to daily means for legibility. Site-level
# R-squared values elsewhere in this workflow remain calculated at 3-hourly
# resolution on the original held-out observations.
permafrost_time_series <- permafrost_representative[
  ,
  .(
    `NEON observations` = mean(Tsoil_obs_C, na.rm = TRUE),
    `Original SIPNET` = mean(Tsoil_original_C, na.rm = TRUE),
    `Optimized tau` = mean(Tsoil_pred_C, na.rm = TRUE)
  ),
  by = .(
    plot_date = as.Date(time)
  )
]


permafrost_delta <- permafrost_sites[
  index == representative_permafrost_index,
  delta_r2
][1L]


permafrost_site_name <- permafrost_sites[
  index == representative_permafrost_index,
  NEON_code
][1L]


permafrost_time_series[
  ,
  panel := sprintf(
    "Permafrost: %s | held-out %s | ΔR² = %.3f",
    permafrost_site_name,
    permafrost_year,
    permafrost_delta
  )
]


representative_time_series <- data.table::rbindlist(
  list(non_time_series, permafrost_time_series),
  use.names = TRUE,
  fill = TRUE
)


representative_time_series <- data.table::melt(
  representative_time_series,
  id.vars = c("plot_date", "panel"),
  measure.vars = c(
    "NEON observations",
    "Original SIPNET",
    "Optimized tau"
  ),
  variable.name = "series",
  value.name = "soil_temperature_C"
)


representative_time_series[
  ,
  series := factor(
    series,
    levels = c(
      "NEON observations",
      "Original SIPNET",
      "Optimized tau"
    )
  )
]


time_series_plot <- ggplot2::ggplot(
  representative_time_series,
  ggplot2::aes(
    x = plot_date,
    y = soil_temperature_C,
    color = series,
    linewidth = series
  )
) +
  ggplot2::geom_line(alpha = 0.9) +
  ggplot2::scale_color_manual(
    values = c(
      "NEON observations" = "#111111",
      "Original SIPNET" = "#CC79A7",
      "Optimized tau" = "#009E73"
    ),
    name = NULL
  ) +
  ggplot2::scale_linewidth_manual(
    values = c(
      "NEON observations" = 0.9,
      "Original SIPNET" = 0.6,
      "Optimized tau" = 0.7
    ),
    guide = "none"
  ) +
  ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x") +
  ggplot2::labs(
    title = "Representative held-out SoilT time series",
    subtitle = paste0(
      "Each site is closest to its model group's median ΔR² | NEON VER",
      target_vertical_position,
      " (~",
      target_depth_cm,
      " cm)"
    ),
    x = "Date",
    y = "Soil temperature (°C)",
    caption = paste0(
      "The permafrost series is shown as daily means for readability; ",
      "its reported validation metrics use the original 3-hour observations."
    )
  ) +
  meeting_theme(14)


figure_manifest[[length(figure_manifest) + 1L]] <- save_figure(
  time_series_plot,
  paste0(
    "10_representative_heldout_tsoil_timeseries_VER",
    target_vertical_position
  ),
  width = 14,
  height = 9
)


# =============================================================================
# Summary outputs
# =============================================================================

figure_manifest_table <- data.table::rbindlist(
  figure_manifest,
  use.names = TRUE,
  fill = TRUE
)


data.table::fwrite(
  figure_manifest_table,
  file.path(configuration$output_dir, "figure_manifest.csv")
)


model_summary <- site_metrics[
  ,
  .(
    n_sites = .N,
    n_sites_with_delta_r2 = sum(is.finite(delta_r2)),
    median_tau_days = stats::median(tau_days, na.rm = TRUE),
    mean_original_r2 = mean(original_r2, na.rm = TRUE),
    mean_improved_r2 = mean(improved_r2, na.rm = TRUE),
    mean_delta_r2 = mean(delta_r2, na.rm = TRUE),
    median_delta_r2 = stats::median(delta_r2, na.rm = TRUE),
    fraction_sites_improved = mean(delta_r2 > 0, na.rm = TRUE),
    mean_delta_rmse_C = mean(delta_rmse_C, na.rm = TRUE)
  ),
  by = model_group
]


data.table::fwrite(
  model_summary,
  file.path(
    configuration$output_dir,
    paste0("model_summary_VER", target_vertical_position, ".csv")
  )
)


run_summary <- c(
  "SOILT / TAU MEETING FIGURE WORKFLOW COMPLETED",
  paste0("Vertical position: ", target_vertical_position),
  paste0("Nominal depth: ", target_depth_cm, " cm"),
  paste0("Non-permafrost sites: ", nrow(non_sites)),
  paste0("Permafrost sites: ", nrow(permafrost_sites)),
  paste0("Figures generated: ", nrow(figure_manifest_table)),
  paste0("Output directory: ", configuration$output_dir),
  "",
  "Primary regression:",
  raw_equation,
  sprintf("R-squared = %.4f", raw_model_summary$r.squared),
  sprintf("Adjusted R-squared = %.4f", raw_model_summary$adj.r.squared)
)


writeLines(
  run_summary,
  file.path(configuration$output_dir, "RUN_SUMMARY.txt")
)


cat(
  paste(run_summary, collapse = "\n"),
  "\n"
)