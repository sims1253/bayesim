#' @keywords internal
#' @importFrom S7 new_class new_property class_character class_function
#' @importFrom S7 class_numeric class_logical
NULL

VALID_CHECKPOINT_FORMATS <- c("rds")

# S7 class definition for SimulationConfig
SimulationConfig <- S7::new_class(
  name = "SimulationConfig",
  properties = list(
    data_grid = S7::new_property(
      class = S7::new_union(S7::class_data.frame, NULL),
      validator = function(value) {
        if (!is.null(value) && nrow(value) < 1) {
          "data_grid must have at least 1 row"
        }
      }
    ),
    fit_grid = S7::new_property(
      class = S7::new_union(S7::class_data.frame, NULL),
      validator = function(value) {
        if (!is.null(value) && nrow(value) < 1) {
          "fit_grid must have at least 1 row"
        }
      }
    ),
    task_grid = S7::new_property(
      class = S7::new_union(S7::class_data.frame, NULL),
      validator = function(value) {
        if (!is.null(value) && nrow(value) < 1) {
          "task_grid must have at least 1 row"
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
    checkpoint_format = S7::new_property(
      class = S7::class_character,
      validator = function(value) {
        if (length(value) != 1 || is.na(value)) {
          "checkpoint_format must be a single character string"
        }
        if (!value %in% VALID_CHECKPOINT_FORMATS) {
          paste0(
            "checkpoint_format must be one of: ",
            paste(VALID_CHECKPOINT_FORMATS, collapse = ", ")
          )
        }
      }
    ),
    chunk_size = S7::new_property(
      class = S7::class_integer,
      validator = function(value) {
        if (length(value) != 1 || is.na(value) || value < 1) {
          "chunk_size must be a positive integer"
        }
      }
    ),
    retain = S7::new_property(
      class = S7::class_list,
      validator = function(value) {
        msg <- validate_resolved_retention_spec(value)
        if (!isTRUE(msg)) {
          msg
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
#' @param task_grid Optional pre-computed task grid. If provided, overrides
#'   data_grid and fit_grid. Must contain either data_spec/fit_spec list-columns
#'   or data_idx/fit_idx index columns.
#' @param data_generator A function with signature `(data_spec, seed, task_ctx) -> data_bundle`.
#'   Generates data for a single replicate given a data specification row.
#' @param fitter An S7 Fitter object that handles model fitting.
#' @param metrics A list of Metric objects.
#' @param n_replicates Positive integer. Number of replicates per data/fit combination.
#' @param seed Integer. Base seed for reproducible random number generation.
#' @param result_path NULL or character path. If provided, results are saved here.
#' @param checkpoint_format Character scalar. Checkpoint storage format.
#'   Currently only `"rds"` is implemented for checkpoint persistence.
#' @param checkpoint_every Positive integer. Save progress every N tasks.
#' @param chunk_size Positive integer. Maximum number of task results to keep
#'   in memory before forcing a checkpoint write. Defaults to `checkpoint_every`.
#' @param max_in_memory Deprecated alias for `chunk_size`.
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
#'   metrics = list(rmse_metric(), bias_metric()),
#'   n_replicates = 100L,
#'   seed = 42L,
#'   checkpoint_format = "rds"
#' )
#' }
simulation_config <- function(
  data_grid = NULL,
  fit_grid = NULL,
  task_grid = NULL,
  data_generator,
  fitter = NULL,
  metrics = NULL,
  n_replicates = 1L,
  seed,
  result_path = NULL,
  checkpoint_format = c("rds"),
  checkpoint_every = 50L,
  chunk_size = NULL,
  max_in_memory = lifecycle::deprecated(),
  retain = c("metrics", "diagnostics"),
  max_errors = Inf
) {
  if (!is.null(task_grid)) {
    if (!is.data.frame(task_grid)) {
      stop(bayesim_config_error("task_grid must be a data.frame"))
    }
    if (nrow(task_grid) < 1) {
      stop(bayesim_config_error("task_grid must have at least 1 row"))
    }

    has_explicit_specs <- all(c("data_spec", "fit_spec") %in% names(task_grid))
    has_index_specs <- all(c("data_idx", "fit_idx") %in% names(task_grid))

    if (!has_explicit_specs && !has_index_specs) {
      stop(bayesim_config_error(
        paste(
          "task_grid must contain either list-columns data_spec/fit_spec or index columns data_idx/fit_idx"
        )
      ))
    }
  }

  if (
    is.null(task_grid) || !all(c("data_spec", "fit_spec") %in% names(task_grid))
  ) {
    if (!is.data.frame(data_grid)) {
      stop(bayesim_config_error("data_grid must be a data.frame"))
    }
    if (nrow(data_grid) < 1) {
      stop(bayesim_config_error("data_grid must have at least 1 row"))
    }

    if (!is.data.frame(fit_grid)) {
      stop(bayesim_config_error("fit_grid must be a data.frame"))
    }
    if (nrow(fit_grid) < 1) {
      stop(bayesim_config_error("fit_grid must have at least 1 row"))
    }
  } else {
    data_grid <- if (is.null(data_grid)) NULL else data_grid
    fit_grid <- if (is.null(fit_grid)) NULL else fit_grid
  }

  if (!is.function(data_generator)) {
    stop(bayesim_config_error("data_generator must be a function"))
  }
  gen_formals <- names(formals(data_generator))
  required_args <- c("data_spec", "seed", "task_ctx")
  if (length(gen_formals) < 3) {
    stop(bayesim_config_error(paste(
      "data_generator must accept at least 3 arguments:",
      "(data_spec, seed, task_ctx)"
    )))
  }

  if (!is.null(fitter)) {
    if (!S7::S7_inherits(fitter)) {
      stop(bayesim_config_error("fitter must be an S7 Fitter object"))
    }
  }

  resolved_metrics <- resolve_metrics(metrics)

  n_replicates <- as.integer(n_replicates)
  if (length(n_replicates) != 1 || is.na(n_replicates) || n_replicates < 1) {
    stop(bayesim_config_error("n_replicates must be a positive integer >= 1"))
  }

  seed <- as.integer(seed)
  if (length(seed) != 1 || is.na(seed)) {
    stop(bayesim_config_error("seed must be a single integer"))
  }

  if (!is.null(result_path)) {
    if (
      !is.character(result_path) ||
        length(result_path) != 1 ||
        is.na(result_path)
    ) {
      stop(bayesim_config_error("result_path must be NULL or a single character string"))
    }
  }

  checkpoint_format <- match.arg(checkpoint_format, VALID_CHECKPOINT_FORMATS)

  checkpoint_every <- as.integer(checkpoint_every)
  if (
    length(checkpoint_every) != 1 ||
      is.na(checkpoint_every) ||
      checkpoint_every < 1
  ) {
    stop(bayesim_config_error("checkpoint_every must be a positive integer >= 1"))
  }

  if (!is.null(chunk_size) && lifecycle::is_present(max_in_memory)) {
    stop(bayesim_config_error("Use either chunk_size or max_in_memory, not both"))
  }

  if (is.null(chunk_size)) {
    chunk_size <- checkpoint_every
  }
  if (lifecycle::is_present(max_in_memory)) {
    lifecycle::deprecate_warn(
      "1.1",
      "simulation_config(max_in_memory)",
      "simulation_config(chunk_size)"
    )
    chunk_size <- max_in_memory
  }

  chunk_size <- as.integer(chunk_size)
  if (
    length(chunk_size) != 1 ||
      is.na(chunk_size) ||
      chunk_size < 1
  ) {
    stop(bayesim_config_error("chunk_size must be a positive integer >= 1"))
  }

  retain <- resolve_retention_spec(retain)

  if (length(max_errors) != 1 || is.na(max_errors)) {
    stop(bayesim_config_error("max_errors must be a single numeric value"))
  }
  if (!is.infinite(max_errors) && max_errors < 0) {
    stop(bayesim_config_error("max_errors must be Inf or a non-negative number"))
  }

  # Create and return S7 object
  SimulationConfig(
    data_grid = data_grid,
    fit_grid = fit_grid,
    task_grid = task_grid,
    data_generator = data_generator,
    fitter = fitter,
    metrics = resolved_metrics,
    n_replicates = n_replicates,
    seed = seed,
    result_path = result_path,
    checkpoint_every = checkpoint_every,
    checkpoint_format = checkpoint_format,
    chunk_size = chunk_size,
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
#' @export
resolve_metrics <- function(metrics) {
  if (is.null(metrics)) {
    return(list())
  }

  if (S7::S7_inherits(metrics, Metric)) {
    return(list(metrics))
  }

  if (is.character(metrics)) {
    stop(bayesim_config_error(
      paste(
        "metrics must be Metric objects, not character names.",
        "Use metric constructors such as list(rmse_metric(), bias_metric())."
      )
    ))
  }

  if (!is.list(metrics)) {
    stop(bayesim_config_error(
      "metrics must be NULL, a Metric object, or a list of Metric objects"
    ))
  }

  for (i in seq_along(metrics)) {
    metric <- metrics[[i]]
    if (!S7::S7_inherits(metric, Metric)) {
      stop(bayesim_config_error(paste0("metrics[[", i, "]] is not an S7 Metric object")))
    }
  }

  metrics
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
  if (!S7::S7_inherits(config, SimulationConfig)) {
    stop(bayesim_config_error("config must be a SimulationConfig object"))
  }

  # Extract properties that define the simulation identity
  # Exclude: result_path, checkpoint_every (runtime settings)
  spec <- list(
    data_grid = config@data_grid,
    fit_grid = config@fit_grid,
    task_grid = config@task_grid,
    data_generator_spec = capture_function_signature(
      config@data_generator
    ),
    fitter_spec = capture_fitter_spec(config@fitter),
    metrics_spec = capture_metrics_spec(config@metrics),
    n_replicates = config@n_replicates,
    seed = config@seed,
    checkpoint_format = config@checkpoint_format,
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
#' @export
capture_function_signature <- function(fn) {
  if (!is.function(fn)) {
    return(NA_character_)
  }

  tryCatch(
    {
      env <- environment(fn)
      env_name <- if (is.null(env)) "" else environmentName(env)
      reference <- NULL

      namespace_name <- NULL
      if (grepl("^namespace:", env_name)) {
        namespace_name <- sub("^namespace:", "", env_name)
      } else if (nzchar(env_name) && env_name %in% loadedNamespaces()) {
        namespace_name <- env_name
      }

      if (!is.null(namespace_name)) {
        ns <- namespace_name
        ns_env <- asNamespace(ns)
        candidates <- ls(ns_env, all.names = TRUE)
        matches <- candidates[vapply(
          candidates,
          function(name) {
            identical(get(name, envir = ns_env), fn)
          },
          logical(1)
        )]

        if (length(matches) == 1) {
          reference <- list(
            package = ns,
            name = matches[[1]],
            version = as.character(getNamespaceVersion(ns))
          )
        }
      }

      args <- formals(fn)
      body_hash <- digest::digest(capture.output(print(body(fn))))

      list(
        rehydratable = !is.null(reference),
        reference = reference,
        environment = env_name,
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
#' @export
capture_fitter_spec <- function(fitter) {
  if (is.null(fitter)) {
    return(NA)
  }

  if (S7::S7_inherits(fitter)) {
    attrs <- attributes(fitter)
    props <- attrs[setdiff(names(attrs), c("class", "S7_class"))]

    return(list(
      class = class(fitter)[1],
      package_version = namespaced_object_version(class(fitter)[1]),
      properties = props,
      rehydratable = grepl("::", class(fitter)[1], fixed = TRUE)
    ))
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
#' @export
capture_metrics_spec <- function(metrics) {
  if (is.null(metrics)) {
    return(list())
  }

  specs <- lapply(metrics, function(m) {
    if (S7::S7_inherits(m)) {
      attrs <- attributes(m)
      props <- attrs[setdiff(names(attrs), c("class", "S7_class"))]

      list(
        class = class(m)[1],
        package_version = namespaced_object_version(class(m)[1]),
        properties = props,
        rehydratable = grepl("::", class(m)[1], fixed = TRUE)
      )
    } else {
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
  if (!S7::S7_inherits(config, SimulationConfig)) {
    stop(bayesim_config_error("config must be a SimulationConfig object"))
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
  S7::S7_inherits(x, SimulationConfig)
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
  if (!S7::S7_inherits(config, SimulationConfig)) {
    stop(bayesim_config_error("config must be a SimulationConfig object"))
  }

  # Check for required components
  if (is.null(config@fitter)) {
    cli::cli_warn("fitter is NULL - simulation may fail during model fitting")
  }

  if (is.null(config@metrics) || length(config@metrics) == 0) {
    cli::cli_warn("No metrics specified - results will be empty")
  }

  if (is.null(config@data_generator)) {
    stop(bayesim_config_error("data_generator cannot be NULL"))
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
#' @note This function performs light validation (is_simulation_config check only).
#'   Full validation is expected at the entry point (run_simulation/resume_simulation)
#'   before this utility function is called.
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
  if (!S7::S7_inherits(config, SimulationConfig)) {
    stop(bayesim_config_error("config must be a SimulationConfig object"))
  }

  if (!is.null(config@task_grid)) {
    return(nrow(config@task_grid))
  }

  nrow(config@data_grid) * nrow(config@fit_grid) * config@n_replicates
}

namespaced_object_version <- function(class_name) {
  if (!is.character(class_name) || length(class_name) != 1) {
    return(NULL)
  }

  parts <- strsplit(class_name, "::", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    return(NULL)
  }

  as.character(getNamespaceVersion(parts[[1]]))
}
