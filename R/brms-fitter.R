#' @title Brms Fitter
#' @description Fitter implementation for brms models. Extends the abstract Fitter class
#'   with brms-specific configuration properties.
#'
#' @param name Character string identifying the fitter (inherited from Fitter)
#' @param supports_predictions Logical indicating if predictions are supported (inherited)
#' @param supports_log_lik Logical indicating if log-likelihood is supported (inherited)
#' @param supports_loo Logical indicating if LOO-CV is supported (inherited)
#' @param supports_epred Logical indicating if posterior expectation
#'   predictions are supported (inherited)
#' @param backend Character string for Stan backend ("cmdstanr" or "rstan")
#' @param chains Integer number of MCMC chains
#' @param iter Integer total iterations per chain
#' @param warmup Integer warmup iterations per chain
#' @param thin Integer thinning interval
#' @param refresh Integer refresh rate for progress output
#' @param silent Integer verbosity level (0, 1, or 2)
#' @param cores Integer number of cores for parallel processing
#' @param precompile Logical; if TRUE (default), the model bank compiles each
#'   distinct model spec once via `brms::brm(chains = 0)` and reuses the prefit
#'   across tasks via `stats::update(recompile = FALSE)`. Set to FALSE to fall
#'   back to a fresh `brms::brm()` per task. When precompiling, specify priors
#'   explicitly: some brms defaults are derived from the template dataset and
#'   would otherwise remain embedded in the reused compiled model.
#' @param allow_default_priors Logical, default FALSE. When FALSE, precompiled
#'   model banks reject model specs without an explicit prior with a fatal
#'   [bayesim_config_error()]: brms derives data-dependent default priors from
#'   the template dataset, and they stay embedded in the compiled model that
#'   the bank reuses for every task (the whole study would silently be fit
#'   with the template's priors). Set TRUE to permit brms data-derived default
#'   priors to be embedded from the template data (rarely what you want; a
#'   notice is emitted once per run). Ignored when `precompile` is FALSE.
#' @param stan_args Named list of Stan/brms arguments passed through to the fit,
#'   e.g. `list(adapt_delta = 0.95, max_treedepth = 12, init = 0.1, threads = 2)`.
#'   NULL (default) uses brms/Stan defaults.
#'
#' @return An S7 BrmsFitter object
#' @export
BrmsFitter <- S7::new_class(
  "BrmsFitter",
  parent = Fitter,
  properties = list(
    # brms implements all optional capabilities declared by Fitter.
    name = S7::new_property(S7::class_character, default = "brms"),
    supports_predictions = S7::new_property(S7::class_logical, default = TRUE),
    supports_log_lik = S7::new_property(S7::class_logical, default = TRUE),
    supports_loo = S7::new_property(S7::class_logical, default = TRUE),
    supports_epred = S7::new_property(S7::class_logical, default = TRUE),
    backend = S7::new_property(S7::class_character, default = "cmdstanr"),
    chains = S7::new_property(S7::class_integer, default = 4L),
    iter = S7::new_property(S7::class_integer, default = 2000L),
    warmup = S7::new_property(S7::class_integer, default = 1000L),
    thin = S7::new_property(S7::class_integer, default = 1L),
    refresh = S7::new_property(S7::class_integer, default = 0L),
    silent = S7::new_property(S7::class_integer, default = 2L),
    cores = S7::new_property(S7::class_integer, default = 1L),
    precompile = S7::new_property(S7::class_logical, default = TRUE),
    allow_default_priors = S7::new_property(S7::class_logical, default = FALSE),
    stan_args = S7::new_property(
      S7::new_union(S7::class_list, NULL),
      default = NULL
    )
  )
)

# Register methods

#' Extract real warmup/sample/total timings from a brmsfit
#'
#' brms stores fit output as a `stanfit` object regardless of backend (cmdstanr
#' output is converted), so [rstan::get_elapsed_time()] works for both backends.
#' Returns per-chain CPU seconds summed across chains. Falls back to
#' `NA_real_` for warmup/sample (and the provided fallback total) when the
#' stanfit timing is unavailable (e.g. a chains = 0 prefit, or an interrupted run).
#'
#' @param fit A brmsfit object.
#' @param fallback_total Numeric length-1; used for `total` when timing can't be
#'   extracted (e.g. the wall-clock elapsed time from a timer).
#'
#' @return Named list: `list(total, warmup, sample)`, all numeric scalars.
#' @keywords internal
extract_brms_timings <- function(fit, fallback_total) {
  tryCatch(
    {
      if (is.null(fit) || is.null(fit$fit)) {
        return(list(
          total = fallback_total,
          warmup = NA_real_,
          sample = NA_real_
        ))
      }
      stanfit <- fit$fit
      if (inherits(stanfit, "CmdStanMCMC") && is.function(stanfit$time)) {
        elapsed <- stanfit$time()
        chains <- elapsed$chains
        if (
          is.data.frame(chains) &&
            all(c("warmup", "sampling") %in% names(chains))
        ) {
          warmup <- sum(chains$warmup)
          sample <- sum(chains$sampling)
          return(list(
            total = warmup + sample,
            warmup = warmup,
            sample = sample
          ))
        }
      }
      if (!inherits(stanfit, "stanfit")) {
        return(list(
          total = fallback_total,
          warmup = NA_real_,
          sample = NA_real_
        ))
      }
      # rstan::get_elapsed_time works on the stanfit object that brms produces
      # for both the cmdstanr and rstan backends. Returns a chains x 2 matrix
      # (warmup, sample) in per-chain seconds. Returns NULL for a chains = 0
      # prefit (no chains were run).
      elapsed <- rstan::get_elapsed_time(stanfit)
      if (is.null(elapsed) || length(elapsed) == 0L) {
        return(list(
          total = fallback_total,
          warmup = NA_real_,
          sample = NA_real_
        ))
      }
      warmup <- sum(elapsed[, "warmup"])
      sample <- sum(elapsed[, "sample"])
      list(total = warmup + sample, warmup = warmup, sample = sample)
    },
    error = function(e) {
      list(total = fallback_total, warmup = NA_real_, sample = NA_real_)
    }
  )
}

#' Build brms control list from stan_args
#'
#' Maps the user-facing `stan_args` fields (`adapt_delta`, `max_treedepth`)
#' into the `control` list that brms/Stan accept, alongside the `init` and
#' `threads` passthroughs. Returns a named list of extra arguments suitable
#' for splicing into `brms::brm()` / `stats::update()` via `do.call()`, or an
#' empty list when `stan_args` is NULL.
#'
#' @param stan_args Named list or NULL.
#'
#' @return A named list (possibly empty) of brms arguments.
#'
#' @keywords internal
build_stan_args_call <- function(stan_args) {
  out <- list()
  if (is.null(stan_args) || length(stan_args) == 0L) {
    return(out)
  }

  # control list: adapt_delta, max_treedepth (and any other control_* keys)
  control_keys <- c("adapt_delta", "max_treedepth")
  control_args <- stan_args[names(stan_args) %in% control_keys]
  # Allow user-supplied explicit control list too.
  if (!is.null(stan_args$control) && is.list(stan_args$control)) {
    control_args <- c(stan_args$control, control_args)
  }
  if (length(control_args) > 0L) {
    out$control <- control_args[!vapply(control_args, is.null, logical(1))]
  }

  if (!is.null(stan_args$init)) {
    out$init <- stan_args$init
  }
  if (!is.null(stan_args$threads)) {
    out$threads <- stan_args$threads
  }
  out
}

#' Run a fresh brms fit (no model bank)
#'
#' Shared helper for the fallback path (`precompile = FALSE` or bank miss).
#'
#' @keywords internal
run_fresh_brms <- function(
  fitter,
  data_bundle,
  fit_spec,
  seed,
  formula,
  family,
  prior,
  stanvars
) {
  stan_call_args <- build_stan_args_call(fitter@stan_args)

  do.call(
    brms::brm,
    c(
      list(
        formula = formula,
        data = data_bundle$train,
        family = family,
        prior = prior,
        stanvars = stanvars,
        backend = fitter@backend,
        chains = fitter@chains,
        iter = fitter@iter,
        warmup = fitter@warmup,
        thin = fitter@thin,
        refresh = fitter@refresh,
        silent = fitter@silent,
        cores = fitter@cores,
        seed = seed
      ),
      stan_call_args
    )
  )
}

#' Update a prefit with new data (model bank path)
#'
#' Runs `stats::update()` with `recompile = FALSE`. BEFORE updating, verifies
#' that the task data + formula + family produce the same Stan *data structure*
#' as the prefit (via `brms::make_standata()`). brms does NOT warn when
#' `recompile = FALSE` is passed against structurally incompatible data — it
#' silently reuses the compiled binary against a wrong model frame — so this
#' explicit structural comparison is the only reliable guard. Raises a fatal
#' [bayesim_internal_error()] on any mismatch (e.g. different predictor count,
#' new factor levels, missing variables).
#'
#' Note: we compare the Stan *data structure* (`make_standata` field names and
#' the design-matrix column count `K`), NOT `make_stancode`, because the latter
#' embeds data-derived default-prior values (e.g. the intercept's `student_t`
#' location = `mean(y)`) that vary across datasets without affecting whether the
#' compiled binary is reusable.
#'
#' @keywords internal
update_prefit <- function(
  prefit,
  fitter,
  data_bundle,
  seed,
  formula,
  family,
  prior,
  stanvars
) {
  stan_call_args <- build_stan_args_call(fitter@stan_args)

  # F6: bank entries are list(prefit, struct_sig). Unpack; the cached struct_sig
  # avoids recomputing brms::make_standata() on the prefit-side per task.
  cached_struct_sig <- NULL
  if (
    is.list(prefit) &&
      !inherits(prefit, "brmsfit") &&
      "prefit" %in% names(prefit)
  ) {
    cached_struct_sig <- prefit$struct_sig
    prefit <- prefit$prefit
  }

  # Structural-mismatch guard: the Stan data STRUCTURE (field names and design
  # matrix dimensions K) must match between the prefit (compiled from template
  # data) and the task data. A mismatch means the compiled binary cannot be
  # reused and brms would silently mis-fit under recompile = FALSE. Done BEFORE
  # update() so a mismatch fails fast. The prefit-side signature is taken from
  # the model bank cache (computed once in build_model_bank()); only the
  # task-side make_standata() runs per task.
  task_struct <- tryCatch(
    brms::make_standata(
      formula = formula,
      data = data_bundle$train,
      family = family
    ),
    error = function(e) NULL
  )
  if (is.null(cached_struct_sig) && !is.null(task_struct)) {
    # Bank did not cache the signature (old/legacy entry): compute it now.
    prefit_struct <- tryCatch(
      brms::make_standata(
        formula = formula,
        data = prefit$data,
        family = family
      ),
      error = function(e) NULL
    )
    cached_struct_sig <- if (!is.null(prefit_struct)) {
      list(
        fields = names(prefit_struct),
        K = prefit_struct$K,
        X_ncol = ncol(prefit_struct$X),
        levels = .data_levels_sig(prefit$data)
      )
    } else {
      NULL
    }
  }
  if (!is.null(cached_struct_sig) && !is.null(task_struct)) {
    # K = number of regression predictors (including intercept); the canonical
    # signal for binary compatibility. Also compare the full field-name set in
    # case the family induces extra data blocks (e.g. shape, theta for ordinal),
    # and the categorical-column levels, whose order silently reorders
    # coefficients under recompile = FALSE.
    task_sig <- list(
      fields = names(task_struct),
      K = task_struct$K,
      X_ncol = ncol(task_struct$X),
      levels = .data_levels_sig(data_bundle$train)
    )
    if (!identical(cached_struct_sig, task_sig)) {
      stop(bayesim_internal_error(paste(
        "Model-bank structural mismatch: the task data has a different Stan",
        "data structure than the template used for precompilation",
        "(predictor count, model-frame shape, or categorical predictor",
        "levels differ). This means the compiled binary cannot be reused.",
        "Ensure all data_grid rows produce data with the same structure,",
        "or set precompile = FALSE."
      )))
    }
  }

  fit_obj <- do.call(
    stats::update,
    c(
      list(
        object = prefit,
        newdata = data_bundle$train,
        seed = seed,
        recompile = FALSE,
        chains = fitter@chains,
        iter = fitter@iter,
        warmup = fitter@warmup,
        thin = fitter@thin,
        refresh = fitter@refresh,
        silent = fitter@silent,
        cores = fitter@cores
      ),
      stan_call_args
    )
  )

  fit_obj
}

#' @export
S7::method(fit_model, BrmsFitter) <- function(
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
      # Resolve model spec from fit_spec (list from fit_grid row). When the
      # grid stores formula/family/prior/stanvars as list-columns, as.list()
      # on a row wraps each in a length-1 list, so unwrap once here.
      unwrap_if_listcol <- function(x) {
        if (length(x) == 1L && is.list(x) && !"data.frame" %in% class(x)) {
          x[[1]]
        } else {
          x
        }
      }
      formula <- unwrap_if_listcol(fit_spec$formula)
      family <- unwrap_if_listcol(fit_spec$family)
      prior <- fit_spec$prior
      stanvars <- fit_spec$stanvars
      # prior/stanvars may also be list-column-wrapped; unwrap if needed, but
      # preserve brmsprior (a data.frame) as-is.
      if (
        length(prior) == 1L && is.list(prior) && !"data.frame" %in% class(prior)
      ) {
        prior <- prior[[1]]
      }
      if (length(stanvars) == 1L && is.list(stanvars)) {
        stanvars <- stanvars[[1]]
      }

      formula <- resolve_formula(formula)
      family <- resolve_family(family)

      # Look up a prefit in the session model bank (set by run_simulation()
      # and pushed to daemons via mirai::everywhere()).
      model_bank <- get_model_bank()
      prefit <- lookup_prefit(
        model_bank = model_bank,
        formula = formula,
        family = family,
        prior = prior,
        stanvars = stanvars,
        backend = fitter@backend
      )

      # Capture warnings emitted during the fit AND during diagnostics
      # extraction (brms emits its Rhat/divergence convergence warning lazily
      # from summary(fit_obj) inside extract_brms_diagnostics, not from brm()/
      # update() themselves). Muffle so they don't surface as test noise. The
      # <<- assigns into the method-local `warnings` accumulator.
      withCallingHandlers(
        {
          fit_obj <- if (!is.null(prefit)) {
            # Model bank path: reuse compiled binary, fail loud on recompile.
            update_prefit(
              prefit = prefit,
              fitter = fitter,
              data_bundle = data_bundle,
              seed = seed,
              formula = formula,
              family = family,
              prior = prior,
              stanvars = stanvars
            )
          } else {
            # Fallback path: fresh compile (precompile = FALSE or bank miss).
            run_fresh_brms(
              fitter = fitter,
              data_bundle = data_bundle,
              fit_spec = fit_spec,
              seed = seed,
              formula = formula,
              family = family,
              prior = prior,
              stanvars = stanvars
            )
          }

          # Extract diagnostics here, inside the handler, because brms' Rhat/
          # divergence convergence warning is emitted lazily by summary(fit_obj)
          # called within extract_brms_diagnostics().
          diag <- tryCatch(
            extract_brms_diagnostics(fit_obj),
            error = function(e) list()
          )
        },
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      timer$stop()

      # Extract draws (draws do not warn; keep outside the handler).
      draws <- tryCatch(
        brms::as_draws_matrix(fit_obj),
        error = function(e) NULL
      )

      new_fit_result(
        success = TRUE,
        fit = fit_obj,
        draws = draws,
        diagnostics = diag,
        timing = extract_brms_timings(
          fit_obj,
          fallback_total = timer$elapsed()
        ),
        warnings = warnings,
        error = NULL,
        data_bundle = data_bundle
      )
    },
    error = function(e) {
      timer$stop()
      # Fatal errors (e.g. model-bank structural mismatch) propagate and stop
      # the whole run; only recoverable errors become a failed task result.
      if (is_fatal_error(e)) {
        stop(e)
      }
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
S7::method(log_lik_matrix, BrmsFitter) <- function(
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
S7::method(predict_epred, BrmsFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL
) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(NULL)
  }
  # F3: expectation predictions (mu, no observation noise) for r2_loo.
  # brms::posterior_epred returns S x N (draws x observations), which is the
  # orientation loo::E_loo expects.
  data <- newdata %||% fit_result$data_bundle$train
  brms::posterior_epred(fit_result$fit, newdata = data)
}

S7::method(loo_fit, BrmsFitter) <- function(fitter, fit_result) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(NULL)
  }

  ll <- brms::log_lik(fit_result$fit)
  # Chain-aware relative efficiency (A4): consistent with build_loo_context()
  # and brms::loo(). ll is S x N (draws x observations).
  r_eff <- relative_eff_from_chains(fitter, fit_result, ll)
  loo_result <- loo::loo(ll, r_eff = r_eff)

  list(
    elpd = loo_result$estimates["elpd_loo", "Estimate"],
    p_loo = loo_result$estimates["p_loo", "Estimate"],
    elpd_se = loo_result$estimates["elpd_loo", "SE"],
    pareto_k = loo::pareto_k_values(loo_result)
  )
}

#' @export
S7::method(fit_diagnostics, BrmsFitter) <- function(fitter, fit_result) {
  if (!fit_result$success || is.null(fit_result$fit)) {
    return(list())
  }

  extract_brms_diagnostics(fit_result$fit)
}

#' Extract brms diagnostics
#'
#' Computes rhat/ESS extrema over **all** parameters (fixed, group-level,
#' distributional, sigma) via `posterior::summarise_draws`, not just the fixed
#' effects from `summary(fit)` (A3). Divergences and max-treedepth hits come
#' from the sampler diagnostics. `lp__` is excluded as it is not a parameter of
#' interest.
#'
#' @param fit brms fit object
#'
#' @return Named list of diagnostics with `rhat_max`, `ess_bulk_min`,
#'   `ess_tail_min`, `divergent`, `max_treedepth`.
#' @keywords internal
extract_brms_diagnostics <- function(fit) {
  # rhat/ESS extrema over all parameters (excluding lp__).
  draw_summary <- tryCatch(
    {
      draws <- posterior::as_draws_array(fit)
      vars <- setdiff(posterior::variables(draws), "lp__")
      posterior::summarise_draws(
        posterior::subset_draws(draws, variable = vars),
        rhat = posterior::rhat,
        ess_bulk = posterior::ess_bulk,
        ess_tail = posterior::ess_tail
      )
    },
    error = function(e) NULL
  )

  if (is.null(draw_summary) || nrow(draw_summary) == 0) {
    rhat_values <- ess_bulk_values <- ess_tail_values <- NA_real_
  } else {
    rhat_values <- draw_summary$rhat
    ess_bulk_values <- draw_summary$ess_bulk
    ess_tail_values <- draw_summary$ess_tail
  }

  # Divergences
  sampler_diag <- brms::nuts_params(fit)
  divergent <- sum(sampler_diag$value[sampler_diag$Parameter == "divergent__"])

  # max_treedepth: read the control setting actually used for this fit via the
  # stored stan_args, with a documented fallback of 10 (A3). The previous
  # `get("control_args", asNamespace("brms"))` reach-in was brittle across
  # brms versions.
  max_treedepth_limit <- tryCatch(
    {
      stan_args <- fit$fit@stan_args[[1]]
      td <- NULL
      if (!is.null(stan_args$control)) {
        td <- stan_args$control$max_treedepth
      }
      if (is.null(td)) {
        td <- 10L
      }
      td
    },
    error = function(e) 10L
  )
  max_treedepth <- sum(
    sampler_diag$value[sampler_diag$Parameter == "treedepth__"] >=
      max_treedepth_limit
  )

  list(
    rhat_max = max(rhat_values, na.rm = TRUE),
    ess_bulk_min = min(ess_bulk_values, na.rm = TRUE),
    ess_tail_min = min(ess_tail_values, na.rm = TRUE),
    divergent = divergent,
    max_treedepth = max_treedepth
  )
}
