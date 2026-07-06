#' @keywords internal
#' @importFrom S7 new_class new_property class_character class_function
#' @importFrom S7 class_numeric class_logical
NULL

VALID_CHECKPOINT_FORMATS <- c("rds")

# I8: valid formats for the optional summary sidecar. "rds" is the default and
# writes nothing extra (the canonical rds checkpoint carries the summary). When
# set to "parquet", the final summary is ALSO written to
# `<result_path>/summary.parquet` for downstream (pandas/arrow/polars) use. The
# rds checkpoint remains the canonical resume artifact either way.
VALID_SUMMARY_FORMATS <- c("rds", "parquet")

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
        # Check function signature has at least 2 arguments
        args <- names(formals(value))
        if (length(args) < 2) {
          "data_generator must accept at least 2 arguments: (data_spec, task_ctx)"
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
      validator = validate_positive_integer(
        message = "n_replicates must be a positive integer"
      )
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
      validator = validate_positive_integer(
        message = "checkpoint_every must be a positive integer"
      )
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
    ),
    daemon_setup = S7::new_property(
      class = S7::new_union(S7::class_function, NULL),
      default = NULL
    ),
    # Workstream I3: optional adaptive stopping policy. NULL = run all tasks.
    # When non-NULL a list with: estimand (character), measure (one of
    # bias/coverage/emp_se/mse/model_se), target_mcse (numeric > 0),
    # min_reps (integer, default 50), check_every (integer, default 50).
    stop_on = S7::new_property(
      class = S7::new_union(S7::class_any, NULL),
      default = NULL
    ),
    # I8: optional parquet sidecar for the summary. Default "rds" writes nothing
    # extra. "parquet" additionally writes `<result_path>/summary.parquet` for
    # downstream consumption. Runtime policy: excluded from the config fingerprint.
    summary_format = S7::new_property(
      class = S7::class_character,
      default = "rds",
      validator = function(value) {
        if (length(value) != 1 || is.na(value)) {
          "summary_format must be a single character string"
        } else if (!value %in% VALID_SUMMARY_FORMATS) {
          paste0(
            "summary_format must be one of: ",
            paste(VALID_SUMMARY_FORMATS, collapse = ", ")
          )
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
#' @param data_generator A function with signature `(data_spec, task_ctx) -> data_bundle`.
#'   Generates data for a single replicate given a data specification row.
#'   `task_ctx$seed` carries the per-task integer seed for backends that need one.
#' @param fitter An S7 Fitter object that handles model fitting.
#' @param metrics A list of Metric objects.
#' @param n_replicates Positive integer. Number of replicates per data/fit combination.
#' @param seed Integer. Base seed for reproducible random number generation.
#' @param result_path NULL or character path. If provided, results are saved here.
#' @param checkpoint_format Character scalar. Checkpoint storage format.
#'   Currently only `"rds"` is implemented for checkpoint persistence. (B4:
#'   excluded from the config fingerprint — it is runtime policy.)
#' @param checkpoint_every Positive integer. Save progress every N tasks. This
#'   single knob also bounds the number of task results held in memory at once
#'   (B4: the former separate `chunk_size` knob was merged into this).
#' @param retain Character vector. What to retain in results. Must be subset of
#'   `c("metrics", "diagnostics", "draws", "predictions", "fit", "data", "warnings")`.
#'   (B4: excluded from the config fingerprint — changing retention must not
#'   invalidate resume.)
#' @param max_errors Numeric. Maximum errors before stopping. Use `Inf` for no
#'   limit. (B4: excluded from the config fingerprint.)
#' @param daemon_setup Optional function run once per mirai daemon (via
#'   `mirai::everywhere()`) before tasks start, e.g. to configure cmdstan
#'   paths or load a model bank. Ignored when no daemons are set. Default NULL.
#' @param stop_on Optional adaptive stopping policy (experimental). `NULL`
#'   (default) runs all tasks. Otherwise a list with elements: `estimand`
#'   (character parameter name), `measure` (one of `"bias"`, `"coverage"`,
#'   `"emp_se"`, `"mse"`, `"model_se"`), `target_mcse` (numeric > 0),
#'   `min_reps` (integer, default 50), `check_every` (integer, default 50).
#'   Once the MCSE of `measure` for `estimand` falls below `target_mcse` AND
#'   at least `min_reps` replicates have completed, remaining pending tasks
#'   are marked `"skipped"` and the run stops. (I3: excluded from the config
#'   fingerprint — it is runtime policy.)
#' @param summary_format Character scalar. Output format for the final summary.
#'   `"rds"` (default) writes nothing extra -- the canonical rds checkpoint
#'   carries the summary and remains the resume artifact. `"parquet"`
#'   additionally writes `<result_path>/summary.parquet` using the suggested
#'   `nanoparquet` package, for downstream consumption (pandas, arrow, polars).
#'   (I8: excluded from the config fingerprint -- runtime policy.)
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
#'   metrics = list(pred_rmse_metric(), pred_bias_metric()),
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
  retain = c("metrics", "diagnostics"),
  max_errors = Inf,
  daemon_setup = NULL,
  stop_on = NULL,
  summary_format = c("rds", "parquet")
) {
  if (!is.null(task_grid)) {
    if (!is.data.frame(task_grid)) {
      cli::cli_abort("task_grid must be a data.frame")
    }
    if (nrow(task_grid) < 1) {
      cli::cli_abort("task_grid must have at least 1 row")
    }

    has_explicit_specs <- all(c("data_spec", "fit_spec") %in% names(task_grid))
    has_index_specs <- all(c("data_idx", "fit_idx") %in% names(task_grid))

    if (!has_explicit_specs && !has_index_specs) {
      cli::cli_abort(
        paste(
          "task_grid must contain either list-columns data_spec/fit_spec or index columns data_idx/fit_idx"
        )
      )
    }
  }

  if (
    is.null(task_grid) || !all(c("data_spec", "fit_spec") %in% names(task_grid))
  ) {
    if (!is.data.frame(data_grid)) {
      cli::cli_abort("data_grid must be a data.frame")
    }
    if (nrow(data_grid) < 1) {
      cli::cli_abort("data_grid must have at least 1 row")
    }

    if (!is.data.frame(fit_grid)) {
      cli::cli_abort("fit_grid must be a data.frame")
    }
    if (nrow(fit_grid) < 1) {
      cli::cli_abort("fit_grid must have at least 1 row")
    }
  } else {
    data_grid <- if (is.null(data_grid)) NULL else data_grid
    fit_grid <- if (is.null(fit_grid)) NULL else fit_grid
  }

  # Validate data_generator
  if (!is.function(data_generator)) {
    cli::cli_abort("data_generator must be a function")
  }
  gen_formals <- names(formals(data_generator))
  required_args <- c("data_spec", "task_ctx")
  if (length(gen_formals) < 2) {
    cli::cli_abort(c(
      "data_generator must accept at least 2 arguments",
      "Required signature: (data_spec, task_ctx)"
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
  if (missing(seed)) {
    cli::cli_abort(c(
      "{.arg seed} is required.",
      i = "bayesim derives every task's RNG stream from this one seed, so an explicit value is required for reproducibility."
    ))
  }
  seed <- as.integer(seed)
  if (length(seed) != 1 || is.na(seed)) {
    cli::cli_abort(c(
      "seed must be a single integer.",
      i = "bayesim derives every task's RNG stream from this one seed."
    ))
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

  checkpoint_format <- match.arg(checkpoint_format, VALID_CHECKPOINT_FORMATS)

  # Validate checkpoint_every
  checkpoint_every <- as.integer(checkpoint_every)
  if (
    length(checkpoint_every) != 1 ||
      is.na(checkpoint_every) ||
      checkpoint_every < 1
  ) {
    cli::cli_abort("checkpoint_every must be a positive integer >= 1")
  }

  retain <- resolve_retention_spec(retain)

  # Validate max_errors
  if (length(max_errors) != 1 || is.na(max_errors)) {
    cli::cli_abort("max_errors must be a single numeric value")
  }
  if (!is.infinite(max_errors) && max_errors < 0) {
    cli::cli_abort("max_errors must be Inf or a non-negative number")
  }

  # I3: validate optional adaptive-stopping policy.
  stop_on <- validate_stop_on(stop_on)

  # I8: resolve summary_format.
  summary_format <- match.arg(summary_format, VALID_SUMMARY_FORMATS)

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
    retain = retain,
    max_errors = max_errors,
    daemon_setup = daemon_setup,
    stop_on = stop_on,
    summary_format = summary_format
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
    return(list())
  }

  if (S7::S7_inherits(metrics, Metric)) {
    return(list(metrics))
  }

  if (is.character(metrics)) {
    cli::cli_abort(
      paste(
        "metrics must be Metric objects, not character names.",
        "Use metric constructors such as list(pred_rmse_metric(), pred_bias_metric())."
      )
    )
  }

  if (!is.list(metrics)) {
    cli::cli_abort(
      "metrics must be NULL, a Metric object, or a list of Metric objects"
    )
  }

  for (i in seq_along(metrics)) {
    metric <- metrics[[i]]
    if (!S7::S7_inherits(metric, Metric)) {
      cli::cli_abort("metrics[[{i}]] is not an S7 Metric object")
    }
  }

  metrics
}

# I3: validate the optional adaptive-stopping policy ----------------------

# Valid performance measures (must match those produced by performance_measures).
VALID_STOP_MEASURES <- c("bias", "coverage", "emp_se", "mse", "model_se")

#' Validate the optional adaptive-stopping policy (I3)
#'
#' `NULL` is valid (no adaptive stopping). Otherwise the input must be a list
#' with: `estimand` (character), `measure` (one of `VALID_STOP_MEASURES`),
#' `target_mcse` (numeric > 0), and optional `min_reps` (integer, default 50)
#' and `check_every` (integer, default 50). Returns a normalized list.
#'
#' @param stop_on NULL or a list.
#' @return NULL or a normalized list.
#' @keywords internal
validate_stop_on <- function(stop_on) {
  if (is.null(stop_on)) {
    return(NULL)
  }
  if (!is.list(stop_on)) {
    stop(bayesim_config_error("stop_on must be NULL or a list"))
  }
  need <- c("estimand", "measure", "target_mcse")
  missing <- setdiff(need, names(stop_on))
  if (length(missing) > 0L) {
    stop(bayesim_config_error(paste(
      "stop_on is missing required element(s):",
      paste(missing, collapse = ", ")
    )))
  }
  estimand <- stop_on$estimand
  if (
    !is.character(estimand) ||
      length(estimand) != 1L ||
      is.na(estimand) ||
      !nzchar(estimand)
  ) {
    stop(bayesim_config_error(
      "stop_on$estimand must be a single non-empty character string"
    ))
  }
  measure <- stop_on$measure
  if (
    !is.character(measure) ||
      length(measure) != 1L ||
      is.na(measure) ||
      !measure %in% VALID_STOP_MEASURES
  ) {
    stop(bayesim_config_error(paste(
      "stop_on$measure must be one of:",
      paste(VALID_STOP_MEASURES, collapse = ", ")
    )))
  }
  target_mcse <- stop_on$target_mcse
  if (
    !is.numeric(target_mcse) ||
      length(target_mcse) != 1L ||
      is.na(target_mcse) ||
      target_mcse <= 0
  ) {
    stop(bayesim_config_error(
      "stop_on$target_mcse must be a single positive numeric value"
    ))
  }
  min_reps <- as.integer(stop_on$min_reps %||% 50L)
  if (length(min_reps) != 1L || is.na(min_reps) || min_reps < 1L) {
    stop(bayesim_config_error(
      "stop_on$min_reps must be a single positive integer"
    ))
  }
  check_every <- as.integer(stop_on$check_every %||% 50L)
  if (length(check_every) != 1L || is.na(check_every) || check_every < 1L) {
    stop(bayesim_config_error(
      "stop_on$check_every must be a single positive integer"
    ))
  }
  list(
    estimand = estimand,
    measure = measure,
    target_mcse = target_mcse,
    min_reps = min_reps,
    check_every = check_every
  )
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
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(...)
#' spec <- as_config_spec(config)
#' # spec can now be serialized or hashed
#' }
as_config_spec <- function(config) {
  if (!S7::S7_inherits(config, SimulationConfig)) {
    cli::cli_abort("config must be a SimulationConfig object")
  }

  # Extract properties that define the simulation identity. B4: exclude runtime
  # policy (result_path, checkpoint_every, checkpoint_format, retain, max_errors)
  # — changing retention or error tolerance must not invalidate resume. I3:
  # stop_on (adaptive stopping) is also runtime policy and excluded here.
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
    seed = config@seed
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
#' The fingerprint excludes runtime policy settings (B4):
#' - `result_path`: Output location doesn't affect simulation identity
#' - `checkpoint_every` / `checkpoint_format`: checkpoint cadence/format is runtime optimization
#' - `retain`: retention policy must not invalidate resume
#' - `max_errors`: error tolerance is runtime policy
#'
#' @param config An S7 SimulationConfig object.
#'
#' @return A character string containing the SHA256 hash of the configuration.
#'
#' @keywords internal
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
    cli::cli_abort("config must be a SimulationConfig object")
  }

  spec <- as_config_spec(config)

  # Convert to JSON for stable serialization, then hash. JSON cannot encode
  # arbitrary R objects (e.g. brms formula/family objects in fit_grid), so fall
  # back to a serialized digest when toJSON fails. The serialized digest is
  # still stable across sessions for the same object structure.
  json_str <- tryCatch(
    jsonlite::toJSON(
      spec,
      auto_unbox = TRUE,
      digits = NA,
      null = "null"
    ),
    error = function(e) NULL
  )

  if (is.null(json_str)) {
    return(digest::digest(spec, algo = "sha256"))
  }

  digest::digest(json_str, algo = "sha256", serialize = FALSE)
}

#' Check if Object is a SimulationConfig
#'
#' @param x An object to check.
#'
#' @return TRUE if x is a SimulationConfig, FALSE otherwise.
#'
#' @keywords internal
is_simulation_config <- function(x) {
  S7::S7_inherits(x, SimulationConfig)
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
#' @keywords internal
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
    cli::cli_abort("config must be a SimulationConfig object")
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
