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
#' @param task A task specification list from the task grid.
#' @param config_spec Plain list config spec (for worker transport)
#' @param fitter S7 Fitter object
#' @param metrics List of Metric objects
#' @param retain Character vector of what to retain
#' @return A bayesim_task_result S3 object. If a recoverable error occurs,
#'   returns a failed task result with error information. Fatal errors are
#'   re-thrown and will stop the simulation.
#'
#' @keywords internal
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
  # C1: run_task_safe is TOTAL — it never throws. Fatal errors are captured
  # into a failed task result carrying `error$fatal = TRUE` and the full
  # condition class chain; the controller (execute_tasks) re-raises a
  # reconstructed condition after collecting the batch. This removes the
  # entire cross-boundary mirai condition-restoration machinery.
  rlang::try_fetch(
    run_task(task, config_spec, fitter, metrics, retain),
    error = function(e) {
      # C1: the handler itself must stay total. On mirai daemons without an
      # installed bayesim (source-loaded controller), namespace resolution of
      # package helpers can fail inside this handler, turning a recoverable
      # task error into a fatal transport error. Fall back to a base-R-only
      # capture that preserves every field the controller consumes.
      tryCatch(
        {
          err_info <- capture_error_info(e)
          # Mark fatal conditions so the controller re-raises them after the
          # batch.
          err_info$fatal <- is_fatal_error(e)
          # Preserve the full condition class chain for reconstruction.
          err_info$condition_class <- class(e)
          new_task_result(
            task_id = task$task_id,
            status = "failed",
            timing = list(total = 0),
            error = err_info
          )
        },
        error = function(handler_error) {
          structure(
            list(
              task_id = task$task_id,
              status = "failed",
              metrics = NULL,
              diagnostics = NULL,
              timing = list(total = 0),
              error = list(
                error_class = paste(class(e), collapse = ", "),
                error_message = tryCatch(
                  paste(conditionMessage(e)),
                  error = function(cond) "unknown error"
                ),
                call = tryCatch(
                  {
                    cl <- conditionCall(e)
                    if (is.null(cl)) NULL else deparse(cl)[1]
                  },
                  error = function(cond) NULL
                ),
                traceback = character(0),
                fatal = length(intersect(
                  class(e),
                  c(
                    "bayesim_config_error",
                    "bayesim_contract_error",
                    "bayesim_checkpoint_error",
                    "bayesim_internal_error"
                  )
                )) >
                  0L,
                condition_class = class(e),
                handler_error = tryCatch(
                  paste(
                    "error handler fell back to base capture:",
                    conditionMessage(handler_error)
                  ),
                  error = function(cond) {
                    "error handler fell back to base capture"
                  }
                )
              ),
              warnings = character(0),
              truth = NULL,
              stop_reason = NULL
            ),
            class = "bayesim_task_result"
          )
        }
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
#' @param task A task specification list from the task grid containing:
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
#' @param fitter S7 Fitter object that implements the fit_model(), predict_fit(),
#'   log_lik_matrix(), and loo_fit() methods
#' @param metrics List of S7 Metric objects to compute
#' @param retain Character vector of what to retain in results. Options:
#'   \itemize{
#'     \item "metrics": Always retained (default)
#'     \item "diagnostics": Fit diagnostics (default)
#'     \item "fit": The raw fit object
#'     \item "draws": Posterior draws matrix
#'   }
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
#' @keywords internal
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
  data_result <- rlang::try_fetch(
    {
      # B4: generator signature is (data_spec, task_ctx); task_ctx$seed carries
      # the integer seed for backends that need one.
      data_bundle <- config_spec$data_generator(
        task$data_spec,
        task_ctx
      )
      validate_data_bundle(data_bundle)
      list(success = TRUE, data_bundle = data_bundle)
    },
    error = function(e) {
      # C1: fatal errors propagate so the controller can stop the run.
      if (is_fatal_error(e)) {
        stop(rlang::cnd_entrace(e))
      }
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
  fit_result <- rlang::try_fetch(
    {
      fit_result <- fit_model(
        fitter,
        data_bundle,
        task$fit_spec,
        task_seed,
        task_ctx
      )
      validate_fit_result_interface(fit_result)
      # The public fitter seam is fit_model() + extract_draws(). A backend may
      # return draws eagerly, but when it returns only its native fit object we
      # obtain the canonical draws here so downstream metrics and retention
      # never need backend-specific knowledge.
      if (
        isTRUE(fit_result$success) &&
          is.null(fit_result$draws) &&
          !is.null(fit_result$fit)
      ) {
        fit_result$draws <- tryCatch(
          extract_draws(fitter, fit_result),
          error = function(e) {
            if (is_fatal_error(e)) {
              stop(e)
            }
            stop(bayesim_contract_error(
              "extract_draws() failed after fit_model(): ",
              conditionMessage(e)
            ))
          }
        )
      }
      if (isTRUE(fit_result$success) && !is.null(fit_result$draws)) {
        validate_fitter_draws(fit_result$draws)
      }
      fit_result
    },
    error = function(e) {
      # C1: fatal errors (config/contract/checkpoint/internal) propagate so the
      # controller can stop the run. Only recoverable fit errors are wrapped.
      if (is_fatal_error(e)) {
        stop(rlang::cnd_entrace(e))
      }
      error_trace <- rlang::trace_back()
      # Normalize to bayesim_fit_error if not already a bayesim error
      if (!is_bayesim_error(e)) {
        e <- bayesim_fit_error(
          message = conditionMessage(e),
          call = conditionCall(e)
        )
      }
      e$trace <- error_trace
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
    task_seed,
    retain = retain
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

  fit_result <- apply_fit_retention(
    fit_result,
    task_retain,
    data_bundle = data_bundle
  )

  task_result <- new_task_result(
    task_id = task$task_id,
    status = "success",
    metrics = metrics_result$metrics,
    diagnostics = fit_result$diagnostics,
    timing = list(total = timer$elapsed()),
    warnings = c(fit_result$warnings, metrics_result$warnings),
    error = NULL,
    truth = data_bundle$true_params
  )

  apply_task_retention(
    task_result,
    fit_result,
    data_bundle,
    task_retain,
    predictions = context$predictions
  )
}

# =============================================================================
# Helper Functions
# =============================================================================

#' Build context for metrics computation
#'
#' Precomputes predictions, log_lik, and loo based on metric needs and retained
#' predictions. Only computes requested context that the fitter supports.
#'
#' @param fit_result A bayesim_fit_result object from a successful fit
#' @param fitter S7 Fitter object
#' @param data_bundle A validated data bundle list
#' @param metrics List of S7 Metric objects
#' @param retain Retention specification. Requesting `"predictions"` computes
#'   predictions even when no metric needs them.
#'
#' @return A named list containing any of:
#'   \itemize{
#'     \item `predictions`: Prediction results from the fitter
#'     \item `log_lik`: Pointwise log-likelihood matrix (S x N)
#'     \item `loo`: LOO-CV results
#'     \item `loo_psis`, `loo_psis_ll`, `loo_epred`: PSIS object, pointwise
#'       log-lik, and posterior-expectation predictions backing the LOO
#'       prediction metrics
#'   }
#'
#' @details
#' The function inspects the `needs` property of each metric to determine
#' what context elements are required. It then checks if the fitter supports
#' each capability and computes them. Any errors during computation result
#' in NULL values for that context element.
#'
#' Evaluation data: `predictions` and `log_lik` are computed on the TEST set
#' when `data_bundle$test` is present, otherwise on the training set. Every
#' built-in metric that consumes them (`pred_*`, `elpd_test`, `r2_test`) compares against the test response, so the predictions must be
#' for the test rows. The LOO context is always built on the training set —
#' leave-one-out is in-sample by construction. `loo_epred` is likewise a
#' training-set matrix; when no metric needs `"loo"` (or the fitter lacks LOO
#' support) it is built directly via [predict_epred()] rather than through
#' the LOO context, so declaring `needs = "epred"` alone still delivers it.
#' The PSIS machinery (`loo_psis`, `loo_psis_ll`) rides along with `loo_epred`:
#' it is computed only when some metric declares `"epred"`; a run whose LOO
#' metrics need only the elpd summary (`needs = "loo"` alone, e.g.
#' `elpd_loo_metric()`) pays for the `loo_fit()` summary alone (#69).
#'
#' @keywords internal
build_metric_context <- function(
  fit_result,
  fitter,
  data_bundle,
  metrics,
  seed = NULL,
  retain = character()
) {
  context <- list()

  # Determine what's needed
  all_needs <- unique(unlist(lapply(metrics, function(m) {
    if (S7::S7_inherits(m)) m@needs else character()
  })))
  retained_options <- if (is.list(retain)) {
    unique(unlist(retain, use.names = FALSE))
  } else {
    retain
  }
  if ("predictions" %in% retained_options) {
    all_needs <- unique(c(all_needs, "predictions"))
  }

  # Warn if metrics need features the fitter doesn't support. Per-task runs
  # would otherwise re-warn thousands of times; .warn_once keeps it to one per
  # run (R1b). The same mismatch is also surfaced once by preflight() before the
  # run, so sequential users see it before any task executes.
  if ("predictions" %in% all_needs && !fitter@supports_predictions) {
    .warn_once(
      "unsupported_capability.predictions",
      "Metric requires predictions but fitter does not support them"
    )
  }
  if ("log_lik" %in% all_needs && !fitter@supports_log_lik) {
    .warn_once(
      "unsupported_capability.log_lik",
      "Metric requires log_lik but fitter does not support it"
    )
  }
  if ("loo" %in% all_needs && !fitter@supports_loo) {
    .warn_once(
      "unsupported_capability.loo",
      "Metric requires loo but fitter does not support it"
    )
  }
  if ("epred" %in% all_needs && !fitter@supports_epred) {
    .warn_once(
      "unsupported_capability.epred",
      "Metric requires epred but fitter does not support it",
      i = "epred-based LOO metrics (rmse_loo, r2_loo) will be NA."
    )
  }

  # Predictions and pointwise log-lik are evaluated on the test set when one
  # exists (all consuming metrics compare against the test response); NULL
  # newdata means the fitter's training data.
  eval_data <- data_bundle$test

  if ("predictions" %in% all_needs && fitter@supports_predictions) {
    context$predictions <- tryCatch(
      {
        predictions <- predict_fit(
          fitter,
          fit_result,
          newdata = eval_data,
          seed = seed
        )
        validate_fitter_predictions(
          predictions,
          n_obs = if (is.null(eval_data)) {
            nrow(data_bundle$train)
          } else {
            nrow(eval_data)
          }
        )
        predictions
      },
      # R1c: a fitter that advertises prediction support but fails to predict is
      # a real anomaly. Surface it once instead of silently degrading every
      # prediction metric to NA with no explanation.
      error = function(e) {
        if (is_fatal_error(e)) {
          stop(e)
        }
        .warn_once(
          "predict_fit_failed",
          "predict_fit() failed; prediction metrics will be NA.",
          i = conditionMessage(e)
        )
        NULL
      }
    )
  }

  if ("log_lik" %in% all_needs && fitter@supports_log_lik) {
    context$log_lik <- tryCatch(
      {
        log_lik <- log_lik_matrix(fitter, fit_result, newdata = eval_data)
        validate_fitter_log_lik(
          log_lik,
          n_obs = if (is.null(eval_data)) {
            nrow(data_bundle$train)
          } else {
            nrow(eval_data)
          }
        )
        log_lik
      },
      error = function(e) {
        if (is_fatal_error(e)) {
          stop(e)
        }
        .warn_once(
          "log_lik_failed",
          "log_lik_matrix() failed; log-likelihood metrics will be NA.",
          i = conditionMessage(e)
        )
        NULL
      }
    )
  }

  loo_ctx <- NULL
  if ("loo" %in% all_needs && fitter@supports_loo) {
    loo_ctx <- tryCatch(
      # #69: only metrics declaring "epred" consume the weighted-prediction
      # machinery (loo_psis/loo_psis_ll/loo_epred); a run with elpd_loo alone
      # pays for the loo_fit() summary only.
      build_loo_context(fitter, fit_result, need_psis = "epred" %in% all_needs),
      error = function(e) {
        .warn_once(
          "loo_build_failed",
          "LOO context could not be built; LOO-based metrics will be NA.",
          i = conditionMessage(e)
        )
        NULL
      }
    )
    if (!is.null(loo_ctx)) {
      # F3: PSIS-based prediction machinery for rmse_loo / r2_loo. Built once
      # here so both metrics share it. May be absent (NULL) if no metric
      # declared "epred" (#69), the fitter does not provide epred/log_lik, or
      # the build failed.
      context$loo <- loo_ctx$loo
      context$loo_psis <- loo_ctx$psis
      context$loo_psis_ll <- loo_ctx$log_lik
      context$loo_epred <- loo_ctx$epred
    }
  }

  # #68: epred does not depend on the LOO machinery, but it used to be built
  # only inside the LOO branch above — a metric declaring needs = "epred"
  # without "loo" (or on a fitter without LOO support) never received
  # context$loo_epred and silently NA-degraded with no warn-once explaining
  # it. Build the matrix directly whenever it is still missing and the LOO
  # context never attempted it:
  # - the LOO context was not built at all (no "loo" need, unsupported LOO,
  #   or a failed build), or
  # - it bailed before epred (train-set log-lik failure).
  # predict_epred() is deliberately not retried when the LOO context already
  # attempted it and failed; that failure is explained by the warn-once
  # below.
  if (
    "epred" %in%
      all_needs &&
      is.null(context[["loo_epred"]]) &&
      (is.null(loo_ctx) || !isTRUE(loo_ctx$epred_attempted))
  ) {
    # Pin an exact NULL "loo" binding when a metric asked for "loo" but the
    # LOO context was not built: $loo would otherwise partial-match the
    # loo_epred matrix below and hand a $-reading metric (elpd_loo) the
    # matrix instead of NULL. context["loo"] <- list(NULL) creates the
    # element; $loo <- NULL would not.
    if ("loo" %in% all_needs && is.null(loo_ctx)) {
      context["loo"] <- list(NULL)
    }
    context$loo_epred <- if (!isTRUE(fitter@supports_epred)) {
      NULL # the unsupported_capability.epred warning above covers this
    } else {
      tryCatch(
        {
          # epred is a training-set quantity (LOO is in-sample by
          # construction); outside a complete LOO context there is no
          # log-lik matrix to align draws against.
          epred <- predict_epred(fitter, fit_result)
          validate_fitter_epred(
            epred,
            n_draws = NULL,
            n_obs = nrow(data_bundle$train)
          )
          epred
        },
        error = function(e) {
          # No is_fatal_error() re-stop: a wrong-shaped predict_epred() return
          # raises a contract error and must degrade through this warn-once
          # exactly as it does when built via the LOO context.
          .warn_once(
            "epred_build_failed",
            "predict_epred() failed; epred-based metrics will be NA.",
            i = conditionMessage(e)
          )
          NULL
        }
      )
    }
  }

  # A fitter that declares epred support but produced no matrix at run
  # time is an anomaly like a failing predict_fit(); the unsupported-
  # capability warning above did not fire, so explain the NA degradation.
  # Only the attempted-and-failed-inside-the-LOO-context case reaches this:
  # every other no-matrix shape was either built (or warned) by the direct
  # branch above.
  if (
    "epred" %in%
      all_needs &&
      isTRUE(fitter@supports_epred) &&
      is.null(context[["loo_epred"]]) &&
      !is.null(loo_ctx) &&
      isTRUE(loo_ctx$epred_attempted)
  ) {
    .warn_once(
      "loo_epred_failed",
      paste0(
        "epred matrix unavailable from the LOO context; ",
        "epred-based LOO metrics (rmse_loo, r2_loo) will be NA."
      )
    )
  }

  context
}

#' Build the LOO context (elpd summary + PSIS object + epred)
#'
#' Constructs the full LOO context for F3's rmse_loo / r2_loo metrics. Computes
#' the elpd/p_loo/pareto_k summary (as the legacy `loo_fit()` did), the PSIS
#' object (for `loo::E_loo()` weighted predictions), the pointwise log-likelihood
#' matrix, and the posterior expectation predictions (epred) — all once, shared
#' across metrics.
#'
#' The PSIS object uses `loo::psis(-ll, r_eff)` with per-observation
#' relative-efficiency factors derived from the fit's chain structure via
#' `posterior::as_draws_df(fit)$.chain` (matches brms' internal `r_eff_log_lik`
#' exactly). Falls back to `r_eff = NULL` (with a captured warning) when chain
#' structure is unavailable, which is mathematically valid but slightly less
#' accurate.
#'
#' epred must be the posterior expectation (mu, no observation noise); for brms
#' this is `brms::posterior_epred`. Only fitters with `supports_epred = TRUE`
#' are asked for it via `predict_epred()`; otherwise epred is NULL and the
#' consuming metrics (r2_loo, rmse_loo) degrade to NA.
#'
#' @param need_psis Logical; whether any metric consumes the
#'   weighted-prediction machinery (`loo_psis`/`loo_psis_ll`/`loo_epred`), i.e.
#'   whether any metric declared the `"epred"` need (#69). When FALSE, only the
#'   `loo_fit()` summary is computed: the train-set log-lik matrix, `r_eff`, the
#'   PSIS object, and epred exist solely to feed that machinery, so a run
#'   configuring elpd_loo alone skips them. The `loo_fit()` summary itself is
#'   independent (fitters compute their own log-lik internally).
#'
#' @return A list with elements `loo`, `psis`, `log_lik`, `epred`, and
#'   `epred_attempted` (logical; whether `predict_epred()` was called), or
#'   NULL on failure. `psis`/`log_lik`/`epred` may be individually NULL if
#'   unavailable; when the train-set log-lik matrix fails the function bails
#'   with `epred_attempted = FALSE` so the caller can still build epred
#'   directly (it does not depend on the log-lik).
#' @keywords internal
build_loo_context <- function(fitter, fit_result, need_psis = FALSE) {
  loo_result <- loo_fit(fitter, fit_result)
  if (is.null(loo_result)) {
    return(NULL)
  }

  # The train-set log-lik matrix, r_eff, the PSIS object, and epred only feed
  # the weighted-prediction machinery; with need_psis = FALSE no metric reads
  # them, so skip the work entirely (#69).
  if (!isTRUE(need_psis)) {
    return(list(
      loo = loo_result,
      psis = NULL,
      log_lik = NULL,
      epred = NULL,
      epred_attempted = FALSE
    ))
  }

  # Pointwise log-likelihood matrix (S x N, draws x observations).
  ll <- tryCatch(log_lik_matrix(fitter, fit_result), error = function(e) NULL)
  if (!is.matrix(ll)) {
    return(list(
      loo = loo_result,
      psis = NULL,
      log_lik = NULL,
      epred = NULL,
      epred_attempted = FALSE
    ))
  }

  # Per-observation relative-efficiency via chain structure.
  r_eff <- tryCatch(
    relative_eff_from_chains(fitter, fit_result, ll),
    error = function(e) NULL
  )

  # PSIS object. If r_eff is NULL (no chain info) loo::psis warns and proceeds
  # without the efficiency correction — mathematically valid, slightly less
  # accurate. Capture the warning so it doesn't surface per-task.
  psis_obj <- if (!is.null(r_eff)) {
    tryCatch(loo::psis(-ll, r_eff = r_eff), error = function(e) NULL)
  } else {
    tryCatch(
      suppressWarnings(loo::psis(-ll)),
      error = function(e) NULL
    )
  }

  # Posterior expectation predictions (mu, no noise). Required for r2_loo;
  # rmse_loo also degrades to NA when epred is absent. Gated on the declared
  # capability: fitters with supports_epred = FALSE skip the call entirely. A
  # wrong-shaped return degrades to NULL here (warn-once at the seam) instead
  # of dying as a generic metric error inside loo::E_loo().
  epred <- if (!isTRUE(fitter@supports_epred)) {
    NULL
  } else {
    tryCatch(
      {
        epred <- predict_epred(fitter, fit_result)
        validate_fitter_epred(epred, n_draws = nrow(ll), n_obs = ncol(ll))
        epred
      },
      error = function(e) NULL
    )
  }

  list(
    loo = loo_result,
    psis = psis_obj,
    log_lik = ll,
    epred = epred,
    epred_attempted = isTRUE(fitter@supports_epred)
  )
}

#' Per-observation relative efficiency from the fit's chain structure
#'
#' Mirrors brms' internal `r_eff_log_lik()` by deriving the chain id per draw
#' from `posterior::as_draws_df(fit)$.chain` and passing it to
#' `loo::relative_eff(exp(ll), chain_id = ...)`. `ll` here is S x N (draws x
#' observations, matching brms' log_lik orientation): `relative_eff` wants the
#' draws dimension as rows so chain_id (length S) matches `nrow(ll)`. Returns
#' one r_eff per observation (length N). Returns NULL when chain information is
#' unavailable (e.g. MockFitter or a fitter whose fit lacks a `.chain`
#' variable), in which case the caller falls back to `r_eff = NULL`.
#' @keywords internal
relative_eff_from_chains <- function(fitter, fit_result, ll) {
  fit <- fit_result$fit
  if (is.null(fit)) {
    return(NULL)
  }
  chain_id <- tryCatch(
    posterior::as_draws_df(fit)$.chain,
    error = function(e) NULL
  )
  # ll is S x N (draws x observations); chain_id has length S (one per draw).
  if (is.null(chain_id) || length(chain_id) != nrow(ll)) {
    return(NULL)
  }
  loo::relative_eff(exp(ll), chain_id = chain_id)
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
  n_metrics <- length(metrics)
  if (n_metrics == 0) {
    return(list(metrics = list(), warnings = character()))
  }

  metric_results <- vector("list", n_metrics)
  warning_acc <- character()

  for (i in seq_len(n_metrics)) {
    metric <- metrics[[i]]
    metric_name <- if (S7::S7_inherits(metric)) metric@name else "unknown"
    is_required <- if (S7::S7_inherits(metric)) {
      isTRUE(metric@required)
    } else {
      FALSE
    }

    metric_result <- tryCatch(
      withr::with_seed(
        derive_metric_seed(task_ctx$seed, metric_name),
        compute_metric(metric, fit_result, data_bundle, context, task_ctx)
      ),
      error = function(e) {
        if (is_required) {
          stop(e)
        }
        list(
          value = NA_real_,
          error_class = class(e)[1],
          error_message = conditionMessage(e)
        )
      }
    )

    if (!is.null(metric_result)) {
      processed <- finalize_metric_values(
        output = metric_result,
        metric_name = metric_name,
        task_ctx = task_ctx,
        result_path = result_path
      )
      metric_results[[i]] <- processed$values
      warning_acc <- c(warning_acc, processed$warnings)
    }
  }

  list(
    metrics = unlist(metric_results, recursive = FALSE),
    warnings = unique(warning_acc)
  )
}

#' Derive a stable, metric-specific RNG seed
#'
#' Isolates stochastic metrics from one another: adding or reordering a metric
#' cannot change another metric's random draws within the same task.
#'
#' @param task_seed Scalar task seed.
#' @param metric_name Character metric name.
#' @return A positive scalar integer seed.
#' @keywords internal
derive_metric_seed <- function(task_seed, metric_name) {
  hash <- digest::digest2int(paste(task_seed, metric_name, sep = ":"))
  as.integer((abs(as.double(hash)) %% (.Machine$integer.max - 1)) + 1)
}

#' Derive a scalar seed from an RNG state vector.
#'
#' Extracts a scalar integer seed from the second element of the L'Ecuyer-CMRG
#' RNG state. Element 2 varies well across parallel streams generated by
#' nextRNGStream(), making it a good source of per-task seed diversity.
#'
#' @param rng_seed Integer vector from create_task_rng_streams()
#'
#' @return Positive integer seed
#'
#' @keywords internal
derive_task_seed <- function(rng_seed) {
  if (is.null(rng_seed) || length(rng_seed) < 2) {
    return(1L)
  }

  as.integer((abs(as.double(rng_seed[[2]])) %% (.Machine$integer.max - 1)) + 1)
}

finalize_metric_values <- function(
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
