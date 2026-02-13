#' @title Brms Fitter
#' @description Fitter implementation for brms models.
#'
#' @export
BrmsFitter <- S7::new_class(
  "BrmsFitter",
  parent = Fitter,
  properties = list(
    # brms-specific properties only (name, supports_* are inherited from Fitter)
    backend = S7::new_property(S7::class_character, default = "rstan"),
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
        error = NULL
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

  # Rhat and ESS
  rhat_values <- summary$fixed[, "Rhat"]
  ess_bulk_values <- summary$fixed[, "Bulk_ESS"]
  ess_tail_values <- summary$fixed[, "Tail_ESS"]

  # Divergences
  sampler_diag <- brms::nuts_params(fit)
  divergent <- sum(sampler_diag$value[sampler_diag$Parameter == "divergent__"])
  max_treedepth <- sum(
    sampler_diag$value[sampler_diag$Parameter == "treedepth__"] >=
      brms::control_args(fit)$max_treedepth
  )

  list(
    rhat_max = max(rhat_values, na.rm = TRUE),
    ess_bulk_min = min(ess_bulk_values, na.rm = TRUE),
    ess_tail_min = min(ess_tail_values, na.rm = TRUE),
    divergent = divergent,
    max_treedepth = max_treedepth
  )
}
