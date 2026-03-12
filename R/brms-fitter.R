#' @title Brms Fitter
#' @description Fitter implementation for brms models.
#'
#' @param name Character string identifying the fitter (inherited from Fitter)
#' @param supports_predictions Logical indicating if predictions are supported (inherited)
#' @param supports_log_lik Logical indicating if log-likelihood is supported (inherited)
#' @param supports_loo Logical indicating if LOO-CV is supported (inherited)
#' @param backend Character string for Stan backend ("cmdstanr" or "rstan")
#' @param chains Integer number of MCMC chains
#' @param iter Integer total iterations per chain
#' @param warmup Integer warmup iterations per chain
#' @param thin Integer thinning interval
#' @param refresh Integer refresh rate for progress output
#' @param silent Integer verbosity level (0, 1, or 2)
#' @param cores Integer number of cores for parallel processing
#'
#' @return An S7 BrmsFitter object
#' @export
BrmsFitter <- S7::new_class(
  "BrmsFitter",
  parent = Fitter,
  properties = list(
    # brms-specific properties only (name, supports_* are inherited from Fitter)
    backend = S7::new_property(S7::class_character, default = "cmdstanr"),
    chains = S7::new_property(S7::class_integer, default = 4L),
    iter = S7::new_property(S7::class_integer, default = 2000L),
    warmup = S7::new_property(S7::class_integer, default = 1000L),
    thin = S7::new_property(S7::class_integer, default = 1L),
    refresh = S7::new_property(S7::class_integer, default = 0L),
    silent = S7::new_property(S7::class_integer, default = 2L),
    cores = S7::new_property(S7::class_integer, default = 1L)
  )
)

# Register methods

#' @export
S7::method(fit, BrmsFitter) <- function(
  fitter,
  data_bundle,
  fit_spec,
  seed,
  task_ctx
) {
  timer <- make_timer()
  timer$start()

  warnings <- character()

  result <- tryCatch(
    {
      # Build brms formula from fit_spec
      formula <- fit_spec$formula
      family <- fit_spec$family

      # Capture warnings
      fit_obj <- withCallingHandlers(
        {
          brms::brm(
            formula = formula,
            data = data_bundle$train,
            family = family,
            backend = fitter@backend,
            chains = fitter@chains,
            iter = fitter@iter,
            warmup = fitter@warmup,
            thin = fitter@thin,
            refresh = fitter@refresh,
            silent = fitter@silent,
            cores = fitter@cores,
            seed = seed
          )
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      timer$stop()

      # Extract draws
      draws <- tryCatch(
        brms::as_draws_matrix(fit_obj),
        error = function(e) NULL
      )

      # Extract diagnostics
      diag <- tryCatch(
        extract_brms_diagnostics(fit_obj),
        error = function(e) list()
      )

      new_fit_result(
        success = TRUE,
        fit = fit_obj,
        draws = draws,
        diagnostics = diag,
        timing = list(
          total = timer$elapsed(),
          warmup = fitter@warmup * (timer$elapsed() / fitter@iter),
          sample = (fitter@iter - fitter@warmup) *
            (timer$elapsed() / fitter@iter)
        ),
        warnings = warnings,
        error = NULL,
        data_bundle = data_bundle
      )
    },
    error = function(e) {
      timer$stop()
      new_fit_result(
        success = FALSE,
        error = e,
        timing = list(total = timer$elapsed(), warmup = 0, sample = 0),
        warnings = warnings
      )
    }
  )

  result
}

#' @export
S7::method(extract_draws, BrmsFitter) <- function(
  fitter,
  fit_result,
  variables = NULL
) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(NULL)
  }

  draws <- brms::as_draws_matrix(fit_result$fit)

  if (!is.null(variables)) {
    draws <- draws[, variables, drop = FALSE]
  }

  draws
}

#' @export
S7::method(predict_fit, BrmsFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL,
  seed = NULL
) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(NULL)
  }

  data <- newdata %||% fit_result$data_bundle$train

  preds <- brms::posterior_predict(
    fit_result$fit,
    newdata = data,
    seed = seed
  )

  list(
    predicted_mean = colMeans(preds),
    predicted_samples = preds,
    predicted_sd = apply(preds, 2, sd)
  )
}

#' @export
S7::method(log_lik, BrmsFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL
) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(NULL)
  }

  data <- newdata %||% fit_result$data_bundle$train

  brms::log_lik(fit_result$fit, newdata = data)
}

#' @export
S7::method(loo, BrmsFitter) <- function(fitter, fit_result) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(NULL)
  }

  ll <- brms::log_lik(fit_result$fit)
  loo_result <- loo::loo(ll)

  list(
    elpd = loo_result$estimates["elpd_loo", "Estimate"],
    p_loo = loo_result$estimates["p_loo", "Estimate"],
    elpd_se = loo_result$estimates["elpd_loo", "SE"],
    pareto_k = loo::pareto_k_values(loo_result)
  )
}

#' @export
S7::method(diagnostics, BrmsFitter) <- function(fitter, fit_result) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(list())
  }

  extract_brms_diagnostics(fit_result$fit)
}

#' Extract brms diagnostics
#'
#' @param fit brms fit object
#'
#' @return Named list of diagnostics
extract_brms_diagnostics <- function(fit) {
  summary <- summary(fit)

  # Rhat and ESS - handle models without fixed effects
  if (!is.null(summary$fixed) && nrow(summary$fixed) > 0) {
    rhat_values <- summary$fixed[, "Rhat"]
    ess_bulk_values <- summary$fixed[, "Bulk_ESS"]
    ess_tail_values <- summary$fixed[, "Tail_ESS"]
  } else {
    rhat_values <- NA_real_
    ess_bulk_values <- NA_real_
    ess_tail_values <- NA_real_
  }

  # Divergences
  sampler_diag <- brms::nuts_params(fit)
  divergent <- sum(sampler_diag$value[sampler_diag$Parameter == "divergent__"])

  # max_treedepth - try to get from fit, default to 10
  max_treedepth <- tryCatch(
    {
      # Try to access control_args from brms namespace
      control_args_fn <- get("control_args", envir = asNamespace("brms"))
      td <- control_args_fn(fit)$max_treedepth
      if (is.null(td)) {
        td <- 10
      } # default
      sum(sampler_diag$value[sampler_diag$Parameter == "treedepth__"] >= td)
    },
    error = function(e) NA_integer_
  )

  list(
    rhat_max = max(rhat_values, na.rm = TRUE),
    ess_bulk_min = min(ess_bulk_values, na.rm = TRUE),
    ess_tail_min = min(ess_tail_values, na.rm = TRUE),
    divergent = divergent,
    max_treedepth = max_treedepth
  )
}
