#' @title Simulation Task Execution
#' @description Functions for executing individual simulation tasks with safe
#'   error handling, metric computation, and retention policies.
#' @name worker
#' @keywords internal
NULL

MAX_INLINE_METRIC_VECTOR_LENGTH <- 50L
MAX_INLINE_METRIC_BYTES <- 64 * 1024

# =============================================================================
# Safe Task Wrapper
# =============================================================================

#' Execute a task with safe error handling
#'
#' Wraps task execution in a tryCatch to ensure recoverable errors are converted
#' to task results. Fatal errors (bayesim_config_error, bayesim_contract_error,
#' etc.) are re-thrown to stop the simulation run.
#'
#' @param task A task specification list (from get_task_spec)
#' @param config_spec Plain list config spec (for worker transport)
#' @param fitter S7 Fitter object
#' @param metrics List of Metric objects
#' @param retain Character vector of what to retain
#'
#' @return A bayesim_task_result S3 object. If a recoverable error occurs,
#'   returns a failed task result with error information. Fatal errors are
#'   re-thrown and will stop the simulation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Run a task safely with error handling
#' result <- run_task_safe(
#'   task = task_spec,
#'   config_spec = config_spec,
#'   fitter = my_fitter,
#'   metrics = list(rmse_metric),
#'   retain = c("metrics", "diagnostics")
#' )
#' }
run_task_safe <- function(
  task,
  config_spec,
  fitter,
  metrics,
  retain = c("metrics", "diagnostics")
) {
  tryCatch(
    run_task(task, config_spec, fitter, metrics, retain),
    error = function(e) {
      # Check if this is a fatal error that should stop the run
      if (is_fatal_error(e)) {
        stop(e) # Re-throw fatal errors
      }
      # Convert recoverable errors to failed task results
      new_task_result(
        task_id = task$task_id,
        status = "failed",
        timing = list(total = 0),
        error = capture_error_info(e)
      )
    }
  )
}

# =============================================================================
# Main Worker Function
# =============================================================================

#' Execute a single simulation task
#'
#' Runs one task: generate data, fit model, compute metrics.
#' All errors are captured and converted to task results.
#'
#' @param task A task specification list (from get_task_spec) containing:
#'   \itemize{
#'     \item `task_id`: Unique task identifier string
#'     \item `data_spec`: Data generation specification
#'     \item `fit_spec`: Model fitting specification
#'     \item `rng_seed`: Precomputed RNG seed for this task
#'     \item `task_ctx`: Task context with indices (data_idx, fit_idx, rep_idx)
#'   }
#' @param config_spec Plain list config spec (for worker transport) containing:
#'   \itemize{
#'     \item `data_generator`: Function to generate data
#'     \item Other configuration parameters
#'   }
#' @param fitter S7 Fitter object that implements the fit(), predict_fit(),
#'   log_lik(), and loo() methods
#' @param metrics List of S7 Metric objects to compute
#' @param retain Character vector of what to retain in results. Options:
#'   \itemize{
#'     \item "metrics": Always retained (default)
#'     \item "diagnostics": Fit diagnostics (default)
#'     \item "fit": The raw fit object
#'     \item "draws": Posterior draws matrix
#'   }
#'
#' @return A bayesim_task_result S3 object containing:
#'   \itemize{
#'     \item `task_id`: Task identifier
#'     \item `status`: "success" or "failed"
#'     \item `metrics`: Named list of computed metrics (if successful)
#'     \item `diagnostics`: Named list of fit diagnostics (if successful)
#'     \item `timing`: List with `total` elapsed time
#'     \item `warnings`: Character vector of warning messages
#'     \item `error`: Error information (if failed)
#'   }
#'
#' @details
#' The task execution follows these steps:
#' \enumerate{
#'   \item Set the RNG state for deterministic reproducibility
#'   \item Generate data using the data_generator function
#'   \item Validate the data bundle
#'   \item Fit the model using the fitter
#'   \item Build context for metric computation (predictions, log_lik, loo)
#'   \item Compute all metrics
#'   \item Apply retention policy to manage memory
#' }
#'
#' Error handling:
#' \itemize{
#'   \item Data generation errors are normalized to bayesim_data_error
#'   \item Fit errors are normalized to bayesim_fit_error
#'   \item Metric errors are handled in compute_all_metrics()
#'   \item Fatal errors (config, contract) propagate and stop the simulation
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a task specification
#' task <- list(
#'   task_id = "d001_f001_r00001",
#'   data_spec = list(n = 100, effect_size = 0.5),
#'   fit_spec = list(prior = "normal(0, 1)"),
#'   rng_seed = create_task_rng_streams(42, 1)[[1]],
#'   task_ctx = list(data_idx = 1, fit_idx = 1, rep_idx = 1)
#' )
#'
#' # Run the task
#' result <- run_task(
#'   task = task,
#'   config_spec = config_spec,
#'   fitter = my_fitter,
#'   metrics = list(rmse_metric, coverage_metric),
#'   retain = c("metrics", "diagnostics")
#' )
#'
#' # Check result
#' if (result$status == "success") {
#'   print(result$metrics)
#' } else {
#' print(result$error$error_message)
#' }
#' }
run_task <- function(
  task,
  config_spec,
  fitter,
  metrics,
  retain = c("metrics", "diagnostics")
) {
  timer <- make_timer()
  timer$start()

  # Set RNG for this task
  set_task_rng(task$rng_seed)
  task_seed <- derive_task_seed(task$rng_seed)
  task_ctx <- c(task$task_ctx, list(seed = task_seed))

  # Step 1: Generate data
  data_result <- tryCatch(
    {
      data_bundle <- config_spec$data_generator(
        task$data_spec,
        task_seed,
        task_ctx
      )
      validate_data_bundle(data_bundle)
      list(success = TRUE, data_bundle = data_bundle)
    },
    error = function(e) {
      # Normalize to bayesim_data_error if not already a bayesim error
      if (!is_bayesim_error(e)) {
        e <- bayesim_data_error(
          message = conditionMessage(e),
          call = conditionCall(e)
        )
      }
      list(success = FALSE, error = capture_error_info(e))
    }
  )

  if (!data_result$success) {
    timer$stop()
    return(new_task_result(
      task_id = task$task_id,
      status = "failed",
      timing = list(total = timer$elapsed()),
      error = data_result$error
    ))
  }

  data_bundle <- data_result$data_bundle

  # Step 2: Fit model
  fit_result <- tryCatch(
    {
      fit(fitter, data_bundle, task$fit_spec, task_seed, task_ctx)
    },
    error = function(e) {
      # Normalize to bayesim_fit_error if not already a bayesim error
      if (!is_bayesim_error(e)) {
        e <- bayesim_fit_error(
          message = conditionMessage(e),
          call = conditionCall(e)
        )
      }
      new_fit_result(
        success = FALSE,
        error = e,
        timing = list(total = 0, warmup = 0, sample = 0)
      )
    }
  )

  if (!fit_result$success) {
    timer$stop()
    return(new_task_result(
      task_id = task$task_id,
      status = "failed",
      timing = list(total = timer$elapsed()),
      error = capture_error_info(fit_result$error)
    ))
  }

  # Step 3: Build context for metrics (predictions, log_lik, loo)
  context <- build_metric_context(
    fit_result,
    fitter,
    data_bundle,
    metrics,
    task_seed
  )

  # Step 4: Compute metrics
  metrics_result <- compute_all_metrics(
    fit_result,
    data_bundle,
    context,
    metrics,
    task_ctx,
    result_path = config_spec$result_path
  )

  timer$stop()

  task_retain <- retention_for_task_result(
    retain,
    "success",
    c(fit_result$warnings, metrics_result$warnings)
  )

  fit_result <- apply_retention(fit_result, data_bundle, task_retain)

  task_result <- new_task_result(
    task_id = task$task_id,
    status = "success",
    metrics = metrics_result$metrics,
    diagnostics = fit_result$diagnostics,
    timing = list(total = timer$elapsed()),
    warnings = c(fit_result$warnings, metrics_result$warnings),
    error = NULL
  )

  apply_task_retention(task_result, fit_result, data_bundle, task_retain)
}

# =============================================================================
# Helper Functions
# =============================================================================

#' Build context for metrics computation
#'
#' Precomputes predictions, log_lik, loo based on metric needs.
#' Only computes what is needed by the metrics, and only if the fitter
#' supports it.
#'
#' @param fit_result A bayesim_fit_result object from a successful fit
#' @param fitter S7 Fitter object
#' @param data_bundle A validated data bundle list
#' @param metrics List of S7 Metric objects
#'
#' @return A named list containing any of:
#'   \itemize{
#'     \item `predictions`: Prediction results from the fitter
#'     \item `log_lik`: Pointwise log-likelihood matrix
#'     \item `loo`: LOO-CV results
#'   }
#'
#' @details
#' The function inspects the `needs` property of each metric to determine
#' what context elements are required. It then checks if the fitter supports
#' each capability and computes them. Any errors during computation result
#' in NULL values for that context element.
#'
#' @keywords internal
build_metric_context <- function(
  fit_result,
  fitter,
  data_bundle,
  metrics,
  seed = NULL
) {
  context <- list()

  # Determine what's needed
  all_needs <- unique(unlist(lapply(metrics, function(m) {
    if (S7::S7_inherits(m)) m@needs else character()
  })))

  # Warn if metrics need features the fitter doesn't support
  if ("predictions" %in% all_needs && !fitter@supports_predictions) {
    cli::cli_warn(
      "Metric requires predictions but fitter does not support them"
    )
  }
  if ("log_lik" %in% all_needs && !fitter@supports_log_lik) {
    cli::cli_warn("Metric requires log_lik but fitter does not support it")
  }
  if ("loo" %in% all_needs && !fitter@supports_loo) {
    cli::cli_warn("Metric requires loo but fitter does not support it")
  }

  if ("predictions" %in% all_needs && fitter@supports_predictions) {
    context$predictions <- tryCatch(
      predict_fit(fitter, fit_result, seed = seed),
      error = function(e) NULL
    )
  }

  if ("log_lik" %in% all_needs && fitter@supports_log_lik) {
    context$log_lik <- tryCatch(
      log_lik(fitter, fit_result),
      error = function(e) NULL
    )
  }

  if ("loo" %in% all_needs && fitter@supports_loo) {
    context$loo <- tryCatch(
      loo(fitter, fit_result),
      error = function(e) NULL
    )
  }

  context
}

#' Compute all metrics for a task
#'
#' Iterates through all metrics and computes their values, handling errors
#' gracefully based on each metric's `required` property.
#'
#' @param fit_result A bayesim_fit_result object from a successful fit
#' @param data_bundle A validated data bundle list
#' @param context A list with precomputed values (from build_metric_context)
#' @param metrics List of S7 Metric objects
#' @param task_ctx Task context with identification information
#'
#' @return A named list with:
#'   \itemize{
#'     \item `metrics`: Named list of computed and flattened metric values
#'     \item `warnings`: Character vector of any warning messages
#'   }
#'
#' @details
#' For each metric:
#' \itemize{
#'   \item If the metric has `required = TRUE` and fails, the error is re-thrown
#'   \item If the metric has `required = FALSE` and fails, NA values are returned
#'     with error information stored in `<metric_name>__error_class` and
#'     `<metric_name>__error_message`
#'   \item All metric outputs are flattened using `flatten_metric_output()`
#' }
#'
#' @keywords internal
compute_all_metrics <- function(
  fit_result,
  data_bundle,
  context,
  metrics,
  task_ctx,
  result_path = NULL
) {
  results <- list()
  warnings <- character()

  for (metric in metrics) {
    metric_name <- if (S7::S7_inherits(metric)) metric@name else "unknown"
    is_required <- if (S7::S7_inherits(metric)) {
      isTRUE(metric@required)
    } else {
      FALSE
    }

    metric_result <- tryCatch(
      {
        compute(metric, fit_result, data_bundle, context, task_ctx)
      },
      error = function(e) {
        if (is_required) {
          stop(e) # Re-throw for required metrics
        }
        # Non-required: return NA values with error info
        list(
          value = NA_real_,
          error_class = class(e)[1],
          error_message = conditionMessage(e)
        )
      }
    )

    # Flatten and add to results
    if (!is.null(metric_result)) {
      processed <- process_metric_output(
        output = metric_result,
        metric_name = metric_name,
        task_ctx = task_ctx,
        result_path = result_path
      )
      flattened <- processed$values
      warnings <- c(warnings, processed$warnings)
      results <- c(results, flattened)
    }
  }

  list(metrics = results, warnings = unique(warnings))
}

#' Apply retention policy to fit result
#'
#' Removes large objects from the fit result based on the retention policy
#' to manage memory usage during simulation runs.
#'
#' @param fit_result A bayesim_fit_result object
#' @param data_bundle A data bundle list (currently unused, for future expansion)
#' @param retain Character vector specifying what to retain. Options:
#'   \itemize{
#'     \item "fit": Keep the raw fit object
#'     \item "draws": Keep the posterior draws matrix
#'     \item "diagnostics": Keep the diagnostics list
#'   }
#'
#' @return The fit_result with specified elements set to NULL based on
#'   the retention policy.
#'
#' @details
#' By default, only metrics and diagnostics are retained. Large objects
#' like the raw fit and draws matrix are removed to minimize memory usage.
#' Users can override this by including "fit" or "draws" in the retain vector.
#'
#' @keywords internal
apply_retention <- function(fit_result, data_bundle, retain) {
  # Remove large objects based on retention policy
  if (!"fit" %in% retain) {
    fit_result$fit <- NULL
  }
  if (!"draws" %in% retain) {
    fit_result$draws <- NULL
  }
  if (!"diagnostics" %in% retain) {
    fit_result$diagnostics <- NULL
  }
  # Remove data_bundle if not explicitly retained
  if (!"data" %in% retain) {
    fit_result$data_bundle <- NULL
  }

  fit_result
}

derive_task_seed <- function(rng_seed) {
  if (is.null(rng_seed) || length(rng_seed) < 2) {
    return(1L)
  }

  as.integer((abs(as.double(rng_seed[[2]])) %% (.Machine$integer.max - 1)) + 1)
}

process_metric_output <- function(
  output,
  metric_name,
  task_ctx,
  result_path = NULL
) {
  validate_metric_output(output, metric_name)

  values <- list()
  warnings <- character()

  for (field_name in names(output)) {
    value <- output[[field_name]]

    if (should_externalize_metric_value(value, result_path)) {
      pointer <- externalize_metric_value(
        value = value,
        metric_name = metric_name,
        field_name = field_name,
        task_ctx = task_ctx,
        result_path = result_path
      )

      prefix <- paste0(metric_name, "__", field_name)
      values[[paste0(prefix, "__externalized")]] <- TRUE
      values[[paste0(prefix, "__artifact_path")]] <- pointer$path
      values[[paste0(prefix, "__artifact_hash")]] <- pointer$hash
      values[[paste0(prefix, "__artifact_size")]] <- pointer$size
      values[[paste0(prefix, "__n_values")]] <- length(value)
      warnings <- c(
        warnings,
        sprintf(
          "Externalized high-cardinality metric '%s__%s' for task '%s'",
          metric_name,
          field_name,
          task_ctx$task_id %||% "unknown"
        )
      )
    } else {
      inline_output <- list(value)
      names(inline_output) <- field_name
      values <- c(values, flatten_metric_output(inline_output, metric_name))
    }
  }

  list(values = values, warnings = unique(warnings))
}

should_externalize_metric_value <- function(value, result_path = NULL) {
  if (is.null(result_path) || length(result_path) == 0) {
    return(FALSE)
  }

  is_named_numeric_vector <- is.double(value) &&
    length(value) > 1 &&
    !is.null(names(value)) &&
    all(names(value) != "" & !is.na(names(value)))

  if (!is_named_numeric_vector) {
    return(FALSE)
  }

  length(value) > MAX_INLINE_METRIC_VECTOR_LENGTH ||
    estimate_size(value) > MAX_INLINE_METRIC_BYTES
}

externalize_metric_value <- function(
  value,
  metric_name,
  field_name,
  task_ctx,
  result_path
) {
  artifacts_dir <- file.path(result_path, "artifacts", "metrics")
  artifact_name <- paste(metric_name, field_name, sep = "__")
  externalize_artifact(value, artifacts_dir, task_ctx$task_id, artifact_name)
}
