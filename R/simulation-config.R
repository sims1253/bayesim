#' @keywords internal
#' @importFrom S7 new_class new_property class_character class_function
#' @importFrom S7 class_numeric class_logical S7_validator
NULL

# Valid retain options
VALID_RETAIN_OPTIONS <- c(
  "metrics",
  "diagnostics",
  "draws",
  "predictions",
  "fit",
  "data",
  "warnings"
)

# S7 class definition for SimulationConfig
SimulationConfig <- S7::new_class(
  name = "SimulationConfig",
  properties = list(
    data_grid = S7::new_property(
      class = S7::class_data.frame,
      validator = function(value) {
        if (nrow(value) < 1) {
          "data_grid must have at least 1 row"
        }
      }
    ),
    fit_grid = S7::new_property(
      class = S7::class_data.frame,
      validator = function(value) {
        if (nrow(value) < 1) {
          "fit_grid must have at least 1 row"
        }
      }
    ),
    data_generator = S7::new_property(
      class = S7::class_function,
      validator = function(value) {
        # Check function signature has at least 3 arguments
        args <- names(formals(value))
        if (length(args) < 3) {
          "data_generator must accept at least 3 arguments: (data_spec, seed, task_ctx)"
        }
      }
    ),
    fitter = S7::new_property(
      class = S7::new_union(S7::class_any, NULL),
      validator = function(value) {
        if (!is.null(value) && !S7::S7_inherits(value)) {
          "fitter must be an S7 object or NULL"
        }
      }
    ),
    metrics = S7::new_property(
      class = S7::new_union(S7::class_list, NULL),
      validator = function(value) {
        if (!is.list(value) && !is.null(value)) {
          "metrics must be a list of Metric objects or NULL"
        }
      }
    ),
    n_replicates = S7::new_property(
      class = S7::class_integer,
      validator = function(value) {
        if (length(value) != 1 || is.na(value) || value < 1) {
          "n_replicates must be a positive integer"
        }
      }
    ),
    seed = S7::new_property(
      class = S7::class_integer,
      validator = function(value) {
        if (length(value) != 1 || is.na(value)) {
          "seed must be a single integer"
        }
      }
    ),
    result_path = S7::new_property(
      class = S7::new_union(S7::class_character, NULL),
      validator = function(value) {
        if (!is.null(value) && (length(value) != 1 || is.na(value))) {
          "result_path must be NULL or a single character string"
        }
      }
    ),
    checkpoint_every = S7::new_property(
      class = S7::class_integer,
      validator = function(value) {
        if (length(value) != 1 || is.na(value) || value < 1) {
          "checkpoint_every must be a positive integer"
        }
      }
    ),
    retain = S7::new_property(
      class = S7::class_character,
      validator = function(value) {
        invalid <- setdiff(value, VALID_RETAIN_OPTIONS)
        if (length(invalid) > 0) {
          paste0(
            "retain contains invalid options: ",
            paste(invalid, collapse = ", "),
            ". Valid options: ",
            paste(VALID_RETAIN_OPTIONS, collapse = ", ")
          )
        }
      }
    ),
    max_errors = S7::new_property(
      class = S7::class_numeric,
      validator = function(value) {
        if (length(value) != 1 || is.na(value)) {
          "max_errors must be a single numeric value"
        } else if (!is.infinite(value) && value < 0) {
          "max_errors must be Inf or a non-negative number"
        }
      }
    )
  )
)

#' Create a Simulation Configuration
#'
#' Creates a SimulationConfig S7 object that fully specifies a simulation study.
#' This configuration defines the data generation grid, fitting grid, metrics,
#' and execution parameters.
#'
#' @param data_grid A data.frame with data generation specifications.
#'   Each row represents a distinct data configuration to simulate.
#' @param fit_grid A data.frame with model fitting specifications.
#'   Each row represents a distinct model configuration to fit.
#' @param data_generator A function with signature `(data_spec, seed, task_ctx) -> data_bundle`.
#'   Generates data for a single replicate given a data specification row.
#' @param fitter An S7 Fitter object that handles model fitting.
#' @param metrics A character vector of metric names or a list of Metric objects.
#'   Character names are resolved via the metric registry.
#' @param n_replicates Positive integer. Number of replicates per data/fit combination.
#' @param seed Integer. Base seed for reproducible random number generation.
#' @param result_path NULL or character path. If provided, results are saved here.
#' @param checkpoint_every Positive integer. Save progress every N tasks.
#' @param retain Character vector. What to retain in results. Must be subset of
#'   `c("metrics", "diagnostics", "draws", "predictions", "fit", "data", "warnings")`.
#' @param max_errors Numeric. Maximum errors before stopping. Use `Inf` for no limit.
#'
#' @return An S7 SimulationConfig object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(
#'   data_grid = data.frame(n = c(100, 500), effect = c(0.5, 1.0)),
#'   fit_grid = data.frame(model = c("baseline", "full")),
#'   data_generator = my_data_gen,
#'   fitter = my_fitter,
#'   metrics = c("rmse", "bias"),
#'   n_replicates = 100L,
#'   seed = 42L
#' )
#' }
simulation_config <- function(
  data_grid,
  fit_grid,
  data_generator,
  fitter = NULL,
  metrics = NULL,
  n_replicates = 1L,
  seed,
  result_path = NULL,
  checkpoint_every = 50L,
  retain = c("metrics", "diagnostics"),
  max_errors = Inf
) {
  # Validate data_grid
  if (!is.data.frame(data_grid)) {
    cli::cli_abort("data_grid must be a data.frame")
  }
  if (nrow(data_grid) < 1) {
    cli::cli_abort("data_grid must have at least 1 row")
  }

  # Validate fit_grid
  if (!is.data.frame(fit_grid)) {
    cli::cli_abort("fit_grid must be a data.frame")
  }
  if (nrow(fit_grid) < 1) {
    cli::cli_abort("fit_grid must have at least 1 row")
  }

  # Validate data_generator
  if (!is.function(data_generator)) {
    cli::cli_abort("data_generator must be a function")
  }
  gen_formals <- names(formals(data_generator))
  required_args <- c("data_spec", "seed", "task_ctx")
  if (length(gen_formals) < 3) {
    cli::cli_abort(c(
      "data_generator must accept at least 3 arguments",
      "Required signature: (data_spec, seed, task_ctx)"
    ))
  }

  # Validate fitter
  if (!is.null(fitter)) {
    # Check if fitter is an S7 object (Fitter class check deferred to when Fitter is defined)
    if (!S7::S7_inherits(fitter)) {
      cli::cli_abort("fitter must be an S7 Fitter object")
    }
  }

  # Resolve metrics
  resolved_metrics <- resolve_metrics(metrics)

  # Validate n_replicates
  n_replicates <- as.integer(n_replicates)
  if (length(n_replicates) != 1 || is.na(n_replicates) || n_replicates < 1) {
    cli::cli_abort("n_replicates must be a positive integer >= 1")
  }

  # Validate seed
  seed <- as.integer(seed)
  if (length(seed) != 1 || is.na(seed)) {
    cli::cli_abort("seed must be a single integer")
  }

  # Validate result_path
  if (!is.null(result_path)) {
    if (
      !is.character(result_path) ||
        length(result_path) != 1 ||
        is.na(result_path)
    ) {
      cli::cli_abort("result_path must be NULL or a single character string")
    }
  }

  # Validate checkpoint_every
  checkpoint_every <- as.integer(checkpoint_every)
  if (
    length(checkpoint_every) != 1 ||
      is.na(checkpoint_every) ||
      checkpoint_every < 1
  ) {
    cli::cli_abort("checkpoint_every must be a positive integer >= 1")
  }

  # Validate retain
  invalid_retain <- setdiff(retain, VALID_RETAIN_OPTIONS)
  if (length(invalid_retain) > 0) {
    cli::cli_abort(c(
      "retain contains invalid options: {invalid_retain}",
      "Valid options: {VALID_RETAIN_OPTIONS}"
    ))
  }
  retain <- match.arg(retain, VALID_RETAIN_OPTIONS, several.ok = TRUE)

  # Validate max_errors
  if (length(max_errors) != 1 || is.na(max_errors)) {
    cli::cli_abort("max_errors must be a single numeric value")
  }
  if (!is.infinite(max_errors) && max_errors < 0) {
    cli::cli_abort("max_errors must be Inf or a non-negative number")
  }

  # Create and return S7 object
  SimulationConfig(
    data_grid = data_grid,
    fit_grid = fit_grid,
    data_generator = data_generator,
    fitter = fitter,
    metrics = resolved_metrics,
    n_replicates = n_replicates,
    seed = seed,
    result_path = result_path,
    checkpoint_every = checkpoint_every,
    retain = retain,
    max_errors = max_errors
  )
}

#' Resolve Metrics to List of Metric Objects
#'
#' Internal helper to convert character metric names or Metric objects
#' into a standardized list of Metric objects.
#'
#' @param metrics Character vector of metric names, list of Metric objects,
#'   or NULL.
#'
#' @return A list of Metric objects, or NULL if input was NULL.
#'
#' @keywords internal
resolve_metrics <- function(metrics) {
  if (is.null(metrics)) {
    return(NULL)
  }

  if (is.character(metrics)) {
    # Look up metric names in registry (to be implemented with metric registry)
    # For now, convert to list of character strings as placeholders
    # This will be updated when the Metric class and registry are available
    resolved <- lapply(metrics, function(name) {
      # Placeholder: return a list with the name
      # When Metric class exists: metric_registry_get(name)
      list(name = name)
    })
    return(resolved)
  }

  if (is.list(metrics)) {
    # Already a list, validate each element
    # When Metric class exists: check each element is Metric
    return(metrics)
  }

  cli::cli_abort("metrics must be a character vector or list of Metric objects")
}

#' Convert SimulationConfig to Plain List for Hashing
#'
#' Converts an S7 SimulationConfig object to a plain list suitable for
#' hashing or serialization. Excludes runtime-specific fields like
#' result_path and checkpoint_every.
#'
#' @param config An S7 SimulationConfig object.
#'
#' @return A named list containing the configuration specification.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(...)
#' spec <- as_config_spec(config)
#' # spec can now be serialized or hashed
#' }
as_config_spec <- function(config) {
  if (!S7::S7_inherits(config) || !inherits(config, "SimulationConfig")) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  # Extract properties that define the simulation identity
  # Exclude: result_path, checkpoint_every (runtime settings)
  spec <- list(
    data_grid = config@data_grid,
    fit_grid = config@fit_grid,
    data_generator_signature = capture_function_signature(
      config@data_generator
    ),
    fitter_spec = capture_fitter_spec(config@fitter),
    metrics_spec = capture_metrics_spec(config@metrics),
    n_replicates = config@n_replicates,
    seed = config@seed,
    retain = config@retain,
    max_errors = config@max_errors
  )

  spec
}

#' Capture Function Signature as String
#'
#' Creates a hashable representation of a function's signature.
#'
#' @param fn A function.
#'
#' @return A character string representing the function signature.
#'
#' @keywords internal
capture_function_signature <- function(fn) {
  if (!is.function(fn)) {
    return(NA_character_)
  }

  tryCatch(
    {
      args <- formals(fn)
      body_hash <- digest::digest(capture.output(print(body(fn))))

      list(
        args = names(args),
        body_hash = body_hash
      )
    },
    error = function(e) {
      NA_character_
    }
  )
}

#' Capture Fitter Specification
#'
#' Creates a hashable representation of a fitter's configuration.
#'
#' @param fitter An S7 Fitter object or NULL.
#'
#' @return A list or NA if NULL.
#'
#' @keywords internal
capture_fitter_spec <- function(fitter) {
  if (is.null(fitter)) {
    return(NA)
  }

  if (S7::S7_inherits(fitter)) {
    # Extract plain properties only - never digest the S7 object directly
    spec <- list(
      class = class(fitter)[1],
      name = tryCatch(fitter@name, error = function(e) NA_character_)
    )

    # Add fitter-specific properties if available
    if (!is.null(fitter@supports_predictions)) {
      spec$supports_predictions <- fitter@supports_predictions
    }
    if (!is.null(fitter@supports_log_lik)) {
      spec$supports_log_lik <- fitter@supports_log_lik
    }
    if (!is.null(fitter@supports_loo)) {
      spec$supports_loo <- fitter@supports_loo
    }

    return(spec)
  }

  NA
}

#' Capture Metrics Specification
#'
#' Creates a hashable representation of metrics configuration.
#'
#' @param metrics A list of Metric objects or NULL.
#'
#' @return A list or NA if NULL.
#'
#' @keywords internal
capture_metrics_spec <- function(metrics) {
  if (is.null(metrics)) {
    return(NA)
  }

  specs <- lapply(metrics, function(m) {
    if (S7::S7_inherits(m)) {
      # Extract plain properties only
      list(
        class = class(m)[1],
        name = tryCatch(m@name, error = function(e) NA_character_),
        needs = tryCatch(as.character(m@needs), error = function(e) {
          character()
        }),
        required = tryCatch(isTRUE(m@required), error = function(e) FALSE)
      )
    } else if (is.list(m) && !is.null(m$name)) {
      # Plain list metric spec
      m$name
    } else {
      # Unknown type - don't digest, return NA
      NA_character_
    }
  })

  specs
}

#' Compute Configuration Fingerprint
#'
#' Computes a cryptographic hash of the simulation configuration.
#' The fingerprint uniquely identifies a simulation configuration
#' for caching and deduplication purposes.
#'
#' The fingerprint excludes runtime-specific settings:
#' - `result_path`: Output location doesn't affect simulation identity
#' - `checkpoint_every`: Checkpoint frequency is runtime optimization
#'
#' @param config An S7 SimulationConfig object.
#'
#' @return A character string containing the SHA256 hash of the configuration.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(
#'   data_grid = data.frame(n = 100),
#'   fit_grid = data.frame(model = "baseline"),
#'   data_generator = my_data_gen,
#'   n_replicates = 10L,
#'   seed = 42L
#' )
#'
#' fingerprint <- compute_config_fingerprint(config)
#' # Use fingerprint for caching or deduplication
#' }
compute_config_fingerprint <- function(config) {
  if (!S7::S7_inherits(config) || !inherits(config, "SimulationConfig")) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  spec <- as_config_spec(config)

  # Convert to JSON for stable serialization, then hash
  json_str <- jsonlite::toJSON(
    spec,
    auto_unbox = TRUE,
    digits = NA,
    null = "null"
  )

  digest::digest(json_str, algo = "sha256", serialize = FALSE)
}

#' Check if Object is a SimulationConfig
#'
#' @param x An object to check.
#'
#' @return TRUE if x is a SimulationConfig, FALSE otherwise.
#'
#' @export
is_simulation_config <- function(x) {
  S7::S7_inherits(x) && inherits(x, "SimulationConfig")
}

#' Validate SimulationConfig Completeness
#'
#' Checks that a SimulationConfig has all required components for running.
#'
#' @param config An S7 SimulationConfig object.
#'
#' @return TRUE if valid, otherwise raises an error.
#'
#' @keywords internal
validate_config_completeness <- function(config) {
  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  # Check for required components
  if (is.null(config@fitter)) {
    cli::cli_warn("fitter is NULL - simulation may fail during model fitting")
  }

  if (is.null(config@metrics) || length(config@metrics) == 0) {
    cli::cli_warn("No metrics specified - results will be empty")
  }

  if (is.null(config@data_generator)) {
    cli::cli_abort("data_generator cannot be NULL")
  }

  TRUE
}

#' Get Total Number of Tasks in Configuration
#'
#' Calculates the total number of simulation tasks based on the configuration.
#' Each task is one (data_spec, fit_spec, replicate) combination.
#'
#' @param config An S7 SimulationConfig object.
#'
#' @return Integer. Total number of tasks.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(
#'   data_grid = data.frame(n = c(100, 500)),  # 2 rows
#'   fit_grid = data.frame(model = c("A", "B")),  # 2 rows
#'   n_replicates = 100L,
#'   ...
#' )
#' get_total_tasks(config)  # Returns 400 (2 * 2 * 100)
#' }
get_total_tasks <- function(config) {
  if (!is_simulation_config(config)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  nrow(config@data_grid) * nrow(config@fit_grid) * config@n_replicates
}
