# CmdStanFitter -----------------------------------------------------------
#
# A fitter for user-supplied Stan programs, for researchers who want their own
# Stan models without brms. See D2. cmdstanr caches compiled binaries by file
# hash natively, so there is no model-bank integration; each daemon
# compiles-or-cache-hits on first use.

# The S7 class. The public constructor (below, same name) auto-derives
# supports_* from the log_lik / epred GQ names.
CmdStanFitter_class <- S7::new_class(
  "CmdStanFitter",
  parent = Fitter,
  properties = list(
    name = S7::new_property(S7::class_character, default = "cmdstan"),
    supports_predictions = S7::new_property(S7::class_logical, default = FALSE),
    supports_log_lik = S7::new_property(S7::class_logical, default = FALSE),
    supports_loo = S7::new_property(S7::class_logical, default = FALSE),
    stan_file = S7::new_property(S7::new_union(S7::class_character, NULL), default = NULL),
    stan_code = S7::new_property(S7::new_union(S7::class_character, NULL), default = NULL),
    stan_data = S7::new_property(S7::class_function),
    log_lik_name = S7::new_property(S7::new_union(S7::class_character, NULL), default = NULL),
    epred_name = S7::new_property(S7::new_union(S7::class_character, NULL), default = NULL),
    chains = S7::new_property(S7::class_integer, default = 4L),
    iter_warmup = S7::new_property(S7::class_integer, default = 1000L),
    iter_sampling = S7::new_property(S7::class_integer, default = 1000L),
    adapt_delta = S7::new_property(S7::new_union(S7::class_numeric, NULL), default = NULL),
    max_treedepth = S7::new_property(S7::new_union(S7::class_integer, NULL), default = NULL),
    parallel_chains = S7::new_property(S7::class_integer, default = 1L),
    init = S7::new_property(S7::class_any, default = NULL),
    extra_args = S7::new_property(S7::new_union(S7::class_list, NULL), default = NULL)
  )
)

#' @title CmdStan Fitter (user-supplied Stan programs)
#' @description
#' Run user-supplied Stan programs via cmdstanr without brms. Compilation is
#' cached by cmdstanr (by file hash); each daemon compiles-or-cache-hits on
#' first use, so there is no model-bank integration (D2).
#'
#' The Stan program may declare generated-quantities blocks for `log_lik` (a
#' vector of pointwise log-likelihoods) and optionally `epred` (the expectation
#' `mu`). Declare their names via the `log_lik` / `epred` arguments.
#'
#' **Newdata prediction is out of scope** (raw Stan has no newdata semantics):
#' `supports_predictions` is FALSE unless `epred` is given, in which case
#' `predict_epred` returns the in-sample GQ matrix and `predict_fit` is
#' unsupported. Test-set metrics require brms or a custom fitter.
#'
#' @param stan_file Path to a `.stan` file, or NULL if `stan_code` is supplied.
#' @param stan_code Character string of Stan code (used when `stan_file` is
#'   NULL). Either `stan_file` or `stan_code` is required.
#' @param stan_data Function `function(data_bundle, fit_spec) -> list` mapping a
#'   bayesim data bundle + fit spec to the Stan `data` list. Required.
#' @param log_lik Name of the generated-quantities log-lik vector in the Stan
#'   program, or NULL if the program has no log_lik GQ. When NULL,
#'   `supports_log_lik`/`supports_loo` are FALSE and elpd metrics degrade.
#' @param epred Optional name of an epred/mu GQ matrix or vector. When supplied,
#'   `supports_predictions` is TRUE and `predict_epred` returns the in-sample GQ
#'   matrix (S x N).
#' @param chains Integer number of MCMC chains.
#' @param iter_warmup Integer warmup iterations per chain.
#' @param iter_sampling Integer sampling iterations per chain.
#' @param adapt_delta Numeric adapt_delta control parameter, or NULL for Stan
#'   default.
#' @param max_treedepth Integer max_treedepth control parameter, or NULL.
#' @param parallel_chains Integer number of chains to run in parallel.
#' @param init Passed through to `$sample()` (NULL, "random", "0", or a list).
#' @param ... Additional named arguments passed through to cmdstanr's `$sample()`.
#'
#' @return An S7 `CmdStanFitter` object.
#' @export
#' @seealso [Fitter], [BrmsFitter], [LinearRegressionFitter]
#' @examples
#' \dontrun{
#' fitter <- CmdStanFitter(
#'   stan_file = "model.stan",
#'   stan_data = function(bundle, spec) list(N = nrow(bundle$train), y = bundle$train$y),
#'   log_lik = "log_lik"
#' )
#' }
CmdStanFitter <- function(
  stan_file = NULL,
  stan_code = NULL,
  stan_data,
  log_lik = NULL,
  epred = NULL,
  chains = 4L,
  iter_warmup = 1000L,
  iter_sampling = 1000L,
  adapt_delta = NULL,
  max_treedepth = NULL,
  parallel_chains = 1L,
  init = NULL,
  ...
) {
  if (is.null(stan_file) && is.null(stan_code)) {
    stop(bayesim_config_error("One of stan_file or stan_code must be supplied"))
  }
  if (!is.function(stan_data)) {
    stop(bayesim_config_error("stan_data must be a function(data_bundle, fit_spec) -> list"))
  }
  dots <- list(...)
  CmdStanFitter_class(
    stan_file = stan_file,
    stan_code = stan_code,
    stan_data = stan_data,
    log_lik_name = log_lik,
    epred_name = epred,
    chains = as.integer(chains),
    iter_warmup = as.integer(iter_warmup),
    iter_sampling = as.integer(iter_sampling),
    adapt_delta = adapt_delta,
    max_treedepth = if (is.null(max_treedepth)) NULL else as.integer(max_treedepth),
    parallel_chains = as.integer(parallel_chains),
    init = init,
    extra_args = if (length(dots)) dots else NULL,
    supports_predictions = !is.null(epred),
    supports_log_lik = !is.null(log_lik),
    supports_loo = !is.null(log_lik)
  )
}

# Internal: drop lp__ and the declared GQ variables from a cmdstan variable
# set. posterior::variables() returns element-wise names ("log_lik[1]", ...),
# so a plain setdiff against the bare GQ name would keep every element; match
# the bare name and its indexed elements.
.cmdstan_keep_vars <- function(vars_all, drop_names) {
  drop_names <- drop_names[!vapply(drop_names, is.null, logical(1))]
  drop_names <- unique(c("lp__", unlist(drop_names)))
  pattern <- paste0(
    "^(", paste(vapply(drop_names, .escape_regex, character(1)), collapse = "|"),
    ")(\\[|$)"
  )
  vars_all[!grepl(pattern, vars_all)]
}

.escape_regex <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

# Internal: compile-or-cache a cmdstan_model. Memoised in a session option
# keyed by file/code hash so each daemon compiles once per distinct stan input.
.cmdstan_model <- function(fitter) {
  cached <- getOption("bayesim.cmdstan_models")
  key <- if (!is.null(fitter@stan_file)) {
    digest::digest(list(file = normalizePath(fitter@stan_file, mustWork = FALSE), code = NULL))
  } else {
    digest::digest(list(file = NULL, code = fitter@stan_code))
  }
  if (!is.null(cached) && !is.null(cached[[key]])) {
    return(cached[[key]])
  }
  mod <- if (!is.null(fitter@stan_file)) {
    cmdstanr::cmdstan_model(stan_file = fitter@stan_file)
  } else {
    cmdstanr::cmdstan_model(stan_code = fitter@stan_code)
  }
  cached[[key]] <- mod
  options(bayesim.cmdstan_models = cached)
  mod
}

# fit_model ----------------------------------------------------------------
S7::method(fit_model, CmdStanFitter_class) <- function(
  fitter, data_bundle, fit_spec, seed, task_ctx
) {
  timer <- make_timer()
  timer$start()
  warnings <- character()

  if (is.null(cmdstanr::cmdstan_version())) {
    stop(bayesim_config_error(
      "cmdstan is not installed; run cmdstanr::install_cmdstan()"
    ))
  }

  # Allow per-task stan_file override via a fit_spec column.
  stan_file <- fit_spec$stan_file %||% fitter@stan_file
  stan_data_list <- fitter@stan_data(data_bundle, fit_spec)
  if (!is.list(stan_data_list)) {
    stop(bayesim_contract_error(
      "stan_data must return a list suitable for cmdstanr"
    ))
  }

  mod <- if (!is.null(stan_file) && !identical(stan_file, fitter@stan_file)) {
    cmdstanr::cmdstan_model(stan_file = stan_file)
  } else {
    .cmdstan_model(fitter)
  }

  sample_args <- list(
    data = stan_data_list,
    seed = seed %||% 0L,
    chains = as.integer(fitter@chains),
    parallel_chains = as.integer(fitter@parallel_chains),
    iter_warmup = as.integer(fitter@iter_warmup),
    iter_sampling = as.integer(fitter@iter_sampling),
    refresh = 0L
  )
  # cmdstanr takes sampler controls as direct $sample() arguments (there is no
  # rstan-style `control` list).
  if (!is.null(fitter@adapt_delta)) sample_args$adapt_delta <- fitter@adapt_delta
  if (!is.null(fitter@max_treedepth)) {
    sample_args$max_treedepth <- as.integer(fitter@max_treedepth)
  }
  if (!is.null(fitter@init)) sample_args$init <- fitter@init
  if (!is.null(fitter@extra_args)) sample_args <- c(sample_args, fitter@extra_args)

  fit <- tryCatch(
    {
      withCallingHandlers(
        do.call(mod$sample, sample_args),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
    },
    error = function(e) {
      stop(bayesim_fit_error(message = conditionMessage(e)))
    }
  )

  timer$stop()

  # Build the S x P draws matrix (exclude lp__ and declared GQ names).
  draws_obj <- fit$draws()
  vars_all <- posterior::variables(draws_obj)
  keep_vars <- .cmdstan_keep_vars(
    vars_all, list(fitter@log_lik_name, fitter@epred_name)
  )
  draws_mat <- if (length(keep_vars)) {
    posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = keep_vars))
  } else {
    posterior::as_draws_matrix(draws_obj)
  }

  fit_obj <- list(
    fitter = "cmdstan",
    data_bundle = data_bundle,
    fit = fit,
    stan_data = stan_data_list,
    log_lik_name = fitter@log_lik_name,
    epred_name = fitter@epred_name,
    n_draws = as.integer(fitter@chains * fitter@iter_sampling),
    seed = seed %||% 0L
  )

  diag <- fit_diagnostics(fitter, list(fit = fit_obj, success = TRUE))

  new_fit_result(
    success = TRUE,
    fit = fit_obj,
    draws = draws_mat,
    diagnostics = diag,
    timing = list(total = timer$elapsed(), warmup = NA_real_, sample = NA_real_),
    warnings = warnings,
    error = NULL,
    data_bundle = data_bundle
  )
}

# extract_draws -----------------------------------------------------------
S7::method(extract_draws, CmdStanFitter_class) <- function(
  fitter, fit_result, variables = NULL
) {
  fit_obj <- fit_result$fit
  cmdstan_fit <- fit_obj$fit
  draws_obj <- cmdstan_fit$draws()
  vars_all <- posterior::variables(draws_obj)
  keep <- .cmdstan_keep_vars(vars_all, list(fitter@log_lik_name, fitter@epred_name))
  if (!is.null(variables)) keep <- intersect(keep, variables)
  if (!length(keep)) return(NULL)
  posterior::as_draws_matrix(posterior::subset_draws(draws_obj, variable = keep))
}

# log_lik_matrix ----------------------------------------------------------
S7::method(log_lik_matrix, CmdStanFitter_class) <- function(
  fitter, fit_result, newdata = NULL
) {
  if (is.null(fitter@log_lik_name)) {
    stop(bayesim_config_error(paste(
      "This CmdStanFitter has no log_lik GQ declared. Provide log_lik = '<name>'",
      "in CmdStanFitter() if a metric needs pointwise log-likelihoods."
    )))
  }
  fit_obj <- fit_result$fit
  draws_obj <- fit_obj$fit$draws()
  ll <- posterior::as_draws_matrix(
    posterior::subset_draws(draws_obj, variable = fitter@log_lik_name)
  )
  if (is.null(dim(ll))) return(NULL)
  ll
}

# predict_epred: in-sample GQ matrix (S x N) --------------------------------
S7::method(predict_epred, CmdStanFitter_class) <- function(fitter, fit_result, newdata = NULL) {
  if (is.null(fitter@epred_name)) return(NULL)
  fit_obj <- fit_result$fit
  draws_obj <- fit_obj$fit$draws()
  ep <- posterior::as_draws_matrix(
    posterior::subset_draws(draws_obj, variable = fitter@epred_name)
  )
  if (is.null(dim(ep))) return(NULL)
  ep
}

# loo_fit ------------------------------------------------------------------
S7::method(loo_fit, CmdStanFitter_class) <- function(fitter, fit_result) {
  if (is.null(fitter@log_lik_name)) {
    return(list(elpd = NA_real_, p_loo = NA_real_, elpd_se = NA_real_, pareto_k = numeric()))
  }
  ll <- tryCatch(log_lik_matrix(fitter, fit_result), error = function(e) NULL)
  if (is.null(ll)) {
    return(list(elpd = NA_real_, p_loo = NA_real_, elpd_se = NA_real_, pareto_k = numeric()))
  }
  loo_result <- tryCatch(suppressWarnings(loo::loo(ll)), error = function(e) NULL)
  if (is.null(loo_result)) {
    return(list(elpd = NA_real_, p_loo = NA_real_, elpd_se = NA_real_, pareto_k = numeric()))
  }
  list(
    elpd = loo_result$estimates["elpd_loo", "Estimate"],
    p_loo = loo_result$estimates["p_loo", "Estimate"],
    elpd_se = loo_result$estimates["elpd_loo", "SE"],
    pareto_k = loo::pareto_k_values(loo_result)
  )
}

# fit_diagnostics ----------------------------------------------------------
S7::method(fit_diagnostics, CmdStanFitter_class) <- function(fitter, fit_result) {
  fit_obj <- fit_result$fit
  cmdstan_fit <- fit_obj$fit

  # Sampler diagnostics (divergences, treedepth, ebfmi).
  sampler_diag <- tryCatch(cmdstan_fit$diagnostic_summary(), error = function(e) NULL)
  divergent <- if (!is.null(sampler_diag)) sum(sampler_diag$divergent %||% 0L) else NA_integer_
  max_treedepth <- if (!is.null(sampler_diag)) sum(sampler_diag$max_treedepth %||% 0L) else NA_integer_

  # rhat/ESS extrema over all parameters (excluding lp__ and GQs).
  rhat_max <- NA_real_; ess_bulk_min <- NA_real_; ess_tail_min <- NA_real_
  draw_summary <- tryCatch(
    {
      draws_obj <- cmdstan_fit$draws()
      vars_all <- posterior::variables(draws_obj)
      keep <- .cmdstan_keep_vars(vars_all, list(fitter@log_lik_name, fitter@epred_name))
      if (length(keep)) {
        posterior::summarise_draws(
          posterior::subset_draws(draws_obj, variable = keep),
          rhat = posterior::rhat,
          ess_bulk = posterior::ess_bulk,
          ess_tail = posterior::ess_tail
        )
      } else NULL
    },
    error = function(e) NULL
  )
  if (!is.null(draw_summary) && nrow(draw_summary)) {
    rhat_max <- max(draw_summary$rhat, na.rm = TRUE)
    ess_bulk_min <- min(draw_summary$ess_bulk, na.rm = TRUE)
    ess_tail_min <- min(draw_summary$ess_tail, na.rm = TRUE)
  }

  list(
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    divergent = divergent,
    max_treedepth = max_treedepth
  )
}

# predict_fit: unsupported for raw Stan (no newdata semantics) -------------
S7::method(predict_fit, CmdStanFitter_class) <- function(
  fitter, fit_result, newdata = NULL, seed = NULL
) {
  stop(bayesim_contract_error(paste(
    "CmdStanFitter does not support posterior-predictive sampling.",
    "Raw Stan has no newdata semantics; use an epred GQ via predict_epred,",
    "or use BrmsFitter / a custom fitter for test-set prediction."
  )))
}
