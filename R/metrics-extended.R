# Extended metric library ------------------------------------------------
#
# Additional Metric subclasses covering the simulation-study metric surface.
# Each follows the metrics-built-in.R pattern: S7 class (@keywords internal),
# exported constructor, then the compute() method. Metrics degrade to NA_real_
# when their required inputs are absent rather than erroring.

# Truth-comparing prediction metrics -------------------------------------

#' @title MAE Metric
#' @description Mean absolute error between predictions and the observed
#'   response on the test set (or training set when no test set).
#' @param name Character string naming the metric. Defaults to "mae".
#' @return A `MaeMetric` object.
#' @keywords internal
MaeMetric <- S7::new_class(
  "MaeMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "mae"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname MaeMetric
#' @description Constructor for MaeMetric.
#' @return A `MaeMetric` object.
#' @export
#' @examples
#' mae_metric()
mae_metric <- function(name = "mae") {
  MaeMetric(name = name, needs = "predictions", required = FALSE)
}

S7::method(compute_metric, MaeMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  if (is.null(context$predictions)) return(list(value = NA_real_))
  test_data <- data_bundle$test %||% data_bundle$train
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean
  list(value = mean(abs(predicted - actual)), n_obs = length(actual))
}

#' @title MSE Metric
#' @description Mean squared error between predictions and the observed
#'   response on the test set (or training set when no test set).
#' @param name Character string naming the metric. Defaults to "mse".
#' @return An `MseMetric` object.
#' @keywords internal
MseMetric <- S7::new_class(
  "MseMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "mse"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname MseMetric
#' @description Constructor for MseMetric.
#' @return An `MseMetric` object.
#' @export
#' @examples
#' mse_metric()
mse_metric <- function(name = "mse") {
  MseMetric(name = name, needs = "predictions", required = FALSE)
}

S7::method(compute_metric, MseMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  if (is.null(context$predictions)) return(list(value = NA_real_))
  test_data <- data_bundle$test %||% data_bundle$train
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean
  list(value = mean((predicted - actual)^2), n_obs = length(actual))
}

#' @title Positivity Probability Metric
#' @description For each parameter in `vars_of_interest`, the posterior
#'   probability that the parameter is positive (fraction of draws > 0).
#' @param name Character string naming the metric. Defaults to "pos_prob".
#' @return A `PosProbMetric` object.
#' @keywords internal
PosProbMetric <- S7::new_class(
  "PosProbMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "pos_prob"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname PosProbMetric
#' @description Constructor for PosProbMetric.
#' @return A `PosProbMetric` object.
#' @export
#' @examples
#' pos_prob_metric()
pos_prob_metric <- function(name = "pos_prob") {
  PosProbMetric(name = name, needs = character(), required = FALSE)
}

S7::method(compute_metric, PosProbMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  if (is.null(fit_result$draws)) return(list(mean = NA_real_))
  draws <- fit_result$draws
  voi <- data_bundle$vars_of_interest %||% colnames(draws)
  mapped <- resolve_draw_columns(voi, colnames(draws))
  if (length(mapped) == 0L) return(list(mean = NA_real_))
  probs <- vapply(names(mapped), function(vname) {
    mean(draws[, mapped[[vname]]] > 0)
  }, numeric(1))
  list(mean = mean(probs), by_param = probs)
}

# Posterior summary metrics ----------------------------------------------

#' @title Posterior Summary Metric
#' @description Posterior mean, median, standard deviation, and quantile-based
#'   credible intervals for each parameter in `vars_of_interest`.
#' @param name Character string naming the metric. Defaults to "posterior_summary".
#' @param prob Credible-interval mass. Defaults to 0.95.
#' @return A `PosteriorSummaryMetric` object.
#' @keywords internal
PosteriorSummaryMetric <- S7::new_class(
  "PosteriorSummaryMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "posterior_summary"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE),
    prob = S7::new_property(S7::class_numeric, default = 0.95)
  )
)

#' @rdname PosteriorSummaryMetric
#' @description Constructor for PosteriorSummaryMetric.
#' @param prob Credible-interval mass. Defaults to 0.95.
#' @return A `PosteriorSummaryMetric` object.
#' @export
#' @examples
#' posterior_summary_metric()
posterior_summary_metric <- function(name = "posterior_summary", prob = 0.95) {
  PosteriorSummaryMetric(name = name, needs = character(), required = FALSE, prob = prob)
}

S7::method(compute_metric, PosteriorSummaryMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  if (is.null(fit_result$draws)) {
    return(list(mean = NA_real_, sd = NA_real_, q_lower = NA_real_, q_upper = NA_real_))
  }
  draws <- fit_result$draws
  voi <- data_bundle$vars_of_interest %||% colnames(draws)
  mapped <- resolve_draw_columns(voi, colnames(draws))
  if (length(mapped) == 0L) {
    return(list(mean = NA_real_, sd = NA_real_, q_lower = NA_real_, q_upper = NA_real_))
  }
  sub <- draws[, mapped, drop = FALSE]
  colnames(sub) <- names(mapped)
  lower_q <- (1 - metric@prob) / 2
  upper_q <- 1 - lower_q
  list(
    mean = colMeans(sub),
    median = apply(sub, 2, stats::median),
    sd = apply(sub, 2, stats::sd),
    q_lower = apply(sub, 2, stats::quantile, probs = lower_q),
    q_upper = apply(sub, 2, stats::quantile, probs = upper_q)
  )
}

# Convergence / sampler diagnostic metrics -------------------------------

#' @title Convergence Metric
#' @description Summarizes convergence diagnostics (max R-hat, min ESS,
#'   divergences) from `fit_result$diagnostics`. Returns NA when absent.
#' @param name Character string naming the metric. Defaults to "convergence".
#' @return A `ConvergenceMetric` object.
#' @keywords internal
ConvergenceMetric <- S7::new_class(
  "ConvergenceMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "convergence"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname ConvergenceMetric
#' @description Constructor for ConvergenceMetric.
#' @return A `ConvergenceMetric` object.
#' @export
#' @examples
#' convergence_metric()
convergence_metric <- function(name = "convergence") {
  ConvergenceMetric(name = name, needs = character(), required = FALSE)
}

S7::method(compute_metric, ConvergenceMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  d <- fit_result$diagnostics
  if (is.null(d)) {
    return(list(rhat_max = NA_real_, ess_bulk_min = NA_real_, ess_tail_min = NA_real_, divergent = NA_integer_))
  }
  list(
    rhat_max = as.numeric(d$rhat_max %||% NA_real_),
    ess_bulk_min = as.numeric(d$ess_bulk_min %||% NA_real_),
    ess_tail_min = as.numeric(d$ess_tail_min %||% NA_real_),
    divergent = as.integer(d$divergent %||% NA_integer_)
  )
}

#' @title Sampler Diagnostics Metric
#' @description Surfaces sampler-level diagnostics (divergences, max treedepth)
#'   from `fit_result$diagnostics`.
#' @param name Character string naming the metric. Defaults to "sampler_diagnostics".
#' @return A `SamplerDiagnosticsMetric` object.
#' @keywords internal
SamplerDiagnosticsMetric <- S7::new_class(
  "SamplerDiagnosticsMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "sampler_diagnostics"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname SamplerDiagnosticsMetric
#' @description Constructor for SamplerDiagnosticsMetric.
#' @return A `SamplerDiagnosticsMetric` object.
#' @export
#' @examples
#' sampler_diagnostics_metric()
sampler_diagnostics_metric <- function(name = "sampler_diagnostics") {
  SamplerDiagnosticsMetric(name = name, needs = character(), required = FALSE)
}

S7::method(compute_metric, SamplerDiagnosticsMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  d <- fit_result$diagnostics
  if (is.null(d)) {
    return(list(divergent = NA_integer_, max_treedepth = NA_integer_))
  }
  list(
    divergent = as.integer(d$divergent %||% NA_integer_),
    max_treedepth = as.integer(d$max_treedepth %||% NA_integer_)
  )
}

# SBC rank metric --------------------------------------------------------

#' @title SBC Rank Metric
#' @description Computes Simulation-Based Calibration ranks: for each parameter
#'   in `vars_of_interest`, counts how many (possibly thinned) posterior draws
#'   are below the data-generating `true_params` value. Under correct
#'   calibration the ranks are uniformly distributed on `0..n_ranks-1`.
#'
#'   Autocorrelation in the posterior draws biases SBC rank uniformity (Talts
#'   et al. 2018, §4.1). By default (`thin = "auto"`) the draws are thinned
#'   toward the minimum bulk-ESS across the ranked variables before ranking:
#'   this keeps ~ESS equally spaced draws, restoring near-independent samples
#'   so the rank distribution is comparable to the standard SBC uniformity
#'   test. `n_ranks` (posterior sample size after thinning + 1 possible ranks)
#'   is reported per variable and is required by the SBC diagnostics
#'   (`plot_rank_ecdf`).
#' @param name Character string naming the metric. Defaults to "rank".
#' @param thin Thinning policy. `"auto"` (default) thins toward the min
#'   bulk-ESS across ranked variables; `FALSE` disables thinning (rank over all
#'   draws); an integer is used directly as the thinning stride.
#' @return A `RankMetric` object.
#' @keywords internal
RankMetric <- S7::new_class(
  "RankMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "rank"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE),
    thin = S7::new_property(
      S7::new_union(S7::class_character, S7::class_logical, S7::class_integer, S7::class_double),
      default = "auto")
  )
)

#' @rdname RankMetric
#' @description Constructor for RankMetric.
#' @param thin Thinning policy: `"auto"` (default), `FALSE`, or an integer stride.
#' @return A `RankMetric` object.
#' @export
#' @examples
#' rank_metric()
#' rank_metric(thin = FALSE)
rank_metric <- function(name = "rank", thin = "auto") {
  RankMetric(name = name, needs = character(), required = FALSE, thin = thin)
}

# Auto-thinning stride toward the min bulk-ESS across ranked variable columns.
# Returns a stride >= 1L (1 = no thinning). Uses posterior::ess_bulk per column;
# if unavailable or degenerate, returns 1L.
.rank_thin_stride <- function(draws, mapped) {
  ess <- tryCatch(
    vapply(names(mapped), function(vname) {
      as.numeric(posterior::ess_bulk(draws[, mapped[[vname]]]))
    }, numeric(1)),
    error = function(e) NA_real_
  )
  if (length(ess) == 0L || all(is.na(ess))) return(1L)
  min_ess <- min(ess, na.rm = TRUE)
  if (!is.finite(min_ess) || min_ess < 1L) return(1L)
  stride <- floor(nrow(draws) / min_ess)
  if (!is.finite(stride) || stride < 1L) return(1L)
  as.integer(stride)
}

S7::method(compute_metric, RankMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  if (is.null(fit_result$draws) || is.null(data_bundle$true_params)) {
    # Keep a (deprecated) `mean` field for back-compat with downstream code
    # that still reads it; it is NA when no ranks are computable.
    return(list(mean = NA_real_))
  }
  draws <- fit_result$draws
  true_params <- data_bundle$true_params
  voi <- data_bundle$vars_of_interest %||% names(true_params)
  mapped <- resolve_draw_columns(voi, colnames(draws))
  # Restrict to vars that also have a truth value.
  mapped <- mapped[names(mapped) %in% names(true_params)]
  if (length(mapped) == 0L) return(list(mean = NA_real_))

  # Determine the thinning stride (F4). auto -> toward min bulk-ESS; integer
  # -> direct stride; FALSE -> no thinning (stride 1).
  thin <- metric@thin
  stride <- if (isFALSE(thin)) {
    1L
  } else if (is.character(thin) && thin == "auto") {
    .rank_thin_stride(draws, mapped)
  } else if (is.numeric(thin) && length(thin) == 1L && !is.na(thin)) {
    max(1L, as.integer(thin))
  } else {
    1L
  }
  thinned <- if (stride > 1L) draws[seq.int(1L, nrow(draws), by = stride), , drop = FALSE] else draws
  n_kept <- nrow(thinned)

  ranks <- vapply(names(mapped), function(vname) {
    sum(thinned[, mapped[[vname]]] < true_params[[vname]])
  }, integer(1))
  # n_ranks per variable: post-thinning sample size + 1 possible ranks (0..S).
  n_ranks <- setNames(rep(n_kept + 1L, length(mapped)), names(mapped))
  # Coerce to double (validate_metric_output requires is.double for named
  # numeric vectors of length > 1) and PRESERVE names (vapply names its output
  # via the character indices, but as.numeric() drops them — re-attach).
  by_param_d <- setNames(as.numeric(ranks), names(mapped))
  n_ranks_d <- setNames(as.numeric(n_ranks), names(mapped))
  list(
    n_draws = as.integer(nrow(draws)),
    stride = as.integer(stride),
    by_param = by_param_d,
    n_ranks = n_ranks_d,
    mean = mean(ranks)
  )
}

# R* metric --------------------------------------------------------------

#' @title R* Diagnostic Metric
#' @description The R* classification-tree diagnostic (Talts et al., 2018).
#'   Requires per-chain draw structure not available in the flattened draws
#'   matrix, so returns NA unless chain information is recovered.
#' @param name Character string naming the metric. Defaults to "rstar".
#' @return An `RstarMetric` object.
#' @keywords internal
RstarMetric <- S7::new_class(
  "RstarMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "rstar"),
    needs = S7::new_property(S7::class_character, default = character()),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname RstarMetric
#' @description Constructor for RstarMetric.
#' @return An `RstarMetric` object.
#' @export
#' @examples
#' rstar_metric()
rstar_metric <- function(name = "rstar") {
  RstarMetric(name = name, needs = character(), required = FALSE)
}

S7::method(compute_metric, RstarMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  list(value = NA_real_)
}

# LOO-based metrics ------------------------------------------------------

#' @title ELPD-LOO Metric
#' @description Expected log-predictive-density via PSIS-LOO from `context$loo`.
#' @param name Character string naming the metric. Defaults to "elpd_loo".
#' @return An `ElpdLooMetric` object.
#' @keywords internal
ElpdLooMetric <- S7::new_class(
  "ElpdLooMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "elpd_loo"),
    needs = S7::new_property(S7::class_character, default = "loo"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname ElpdLooMetric
#' @description Constructor for ElpdLooMetric.
#' @return An `ElpdLooMetric` object.
#' @export
#' @examples
#' elpd_loo_metric()
elpd_loo_metric <- function(name = "elpd_loo") {
  ElpdLooMetric(name = name, needs = "loo", required = FALSE)
}

S7::method(compute_metric, ElpdLooMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  loo <- context$loo
  if (is.null(loo)) {
    return(list(elpd = NA_real_, p_loo = NA_real_, se = NA_real_, pareto_k_max = NA_real_))
  }
  list(
    elpd = as.numeric(loo$elpd %||% NA_real_),
    p_loo = as.numeric(loo$p_loo %||% NA_real_),
    se = as.numeric(loo$elpd_se %||% NA_real_),
    pareto_k_max = as.numeric(max(loo$pareto_k, na.rm = TRUE) %||% NA_real_)
  )
}

#' @title RMSE-LOO Metric
#' @description RMSE estimated via PSIS-LOO: the LOO-weighted posterior-mean
#'   prediction is computed with [loo::E_loo()] and compared to the observed
#'   response (`sqrt(mean((y - yloo)^2))`). This matches the construction
#'   underlying brms' `loo_predict(type = "mean")`. The max Pareto-k is
#'   returned as a diagnostic. Falls back to NA when the PSIS object or
#'   observed response is unavailable.
#' @param name Character string naming the metric. Defaults to "rmse_loo".
#' @return A `RmseLooMetric` object.
#' @keywords internal
RmseLooMetric <- S7::new_class(
  "RmseLooMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "rmse_loo"),
    needs = S7::new_property(S7::class_character, default = "loo"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname RmseLooMetric
#' @description Constructor for RmseLooMetric.
#' @return A `RmseLooMetric` object.
#' @export
#' @examples
#' rmse_loo_metric()
rmse_loo_metric <- function(name = "rmse_loo") {
  RmseLooMetric(name = name, needs = "loo", required = FALSE)
}

S7::method(compute_metric, RmseLooMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  loo <- context$loo
  psis_obj <- context$loo_psis
  ll <- context$loo_psis_ll
  ppred <- context$loo_epred
  # F3: PSIS-based LOO-RMSE. Requires the PSIS object + pointwise log-lik,
  # plus a posterior-prediction matrix to weight. We prefer epred (mu) for
  # consistency with brms loo_predict; posterior_predict (with noise) is also
  # a valid mean-type LOO prediction. The precomputed context ships epred; if
  # absent, we cannot compute the LOO mean prediction and degrade to NA.
  if (is.null(loo) || is.null(psis_obj) || is.null(ll) || is.null(ppred)) {
    return(list(value = NA_real_, elpd = NA_real_, pareto_k_max = NA_real_))
  }
  test_data <- data_bundle$test %||% data_bundle$train
  y <- test_data[[data_bundle$response]]
  if (is.null(y)) {
    return(list(value = NA_real_, elpd = NA_real_, pareto_k_max = NA_real_))
  }
  eloo <- loo::E_loo(ppred, psis_obj, log_ratios = -ll, type = "mean")
  yloo <- eloo$value
  pareto_k_max <- suppressWarnings(
    as.numeric(max(eloo$pareto_k, na.rm = TRUE))
  )
  if (!is.numeric(pareto_k_max) || length(pareto_k_max) == 0L || is.infinite(pareto_k_max)) {
    pareto_k_max <- NA_real_
  }
  list(
    elpd = as.numeric(loo$elpd %||% NA_real_),
    value = sqrt(mean((y - yloo)^2)),
    pareto_k_max = pareto_k_max
  )
}

#' @title R-squared LOO Metric
#' @description Variance-explained estimated via PSIS-LOO following brms'
#'   `loo_R2()`: `1 - var_loo(y - yloo) / var_loo(y)`, where `yloo` is the
#'   LOO-weighted posterior expectation ([loo::E_loo()] mean of `posterior_epred`)
#'   and the variances use the same weighted-expecation construction as brms
#'   (Gelman, Goodrich, Gabry & Vehtari 2018). Falls back to NA when epred or
#'   the PSIS object is unavailable.
#' @param name Character string naming the metric. Defaults to "r2_loo".
#' @return An `R2LooMetric` object.
#' @keywords internal
R2LooMetric <- S7::new_class(
  "R2LooMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "r2_loo"),
    needs = S7::new_property(S7::class_character, default = "loo"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname R2LooMetric
#' @description Constructor for R2LooMetric.
#' @return An `R2LooMetric` object.
#' @export
#' @examples
#' r2_loo_metric()
r2_loo_metric <- function(name = "r2_loo") {
  R2LooMetric(name = name, needs = "loo", required = FALSE)
}

S7::method(compute_metric, R2LooMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  loo <- context$loo
  psis_obj <- context$loo_psis
  ll <- context$loo_psis_ll
  epred <- context$loo_epred
  # F3: PSIS-based LOO-R2 (brms loo_R2 construction). Requires epred (mu, no
  # noise), the PSIS object, and the pointwise log-lik. r2_loo MUST use the
  # expectation, not posterior_predict noise draws — see PLAN.md F3.
  if (is.null(loo) || is.null(psis_obj) || is.null(ll) || is.null(epred)) {
    return(list(value = NA_real_, elpd = NA_real_))
  }
  test_data <- data_bundle$test %||% data_bundle$train
  y <- test_data[[data_bundle$response]]
  if (is.null(y)) {
    return(list(value = NA_real_, elpd = NA_real_))
  }
  # E_loo mean of the expectation draws = the LOO point prediction.
  yloo <- loo::E_loo(epred, psis_obj, log_ratios = -ll, type = "mean")$value
  err_loo <- yloo - y
  S <- nrow(epred)
  N <- ncol(epred)
  # brms .loo_R2 variance construction: expectation over Exp(1)-weighted draws.
  exp_draws <- matrix(rexp(S * N, rate = 1), nrow = S, ncol = N)
  weights <- exp_draws / rowSums(exp_draws)
  var_y <- (N / (N - 1)) * (
    rowSums(sweep(weights, 2, y^2, FUN = "*")) -
      rowSums(sweep(weights, 2, y, FUN = "*"))^2
  )
  var_err_loo <- (N / (N - 1)) * (
    rowSums(sweep(weights, 2, err_loo^2, FUN = "*")) -
      rowSums(sweep(weights, 2, err_loo, FUN = "*"))^2
  )
  r2 <- 1 - var_err_loo / var_y
  r2[r2 < -1] <- -1
  r2[r2 > 1] <- 1
  list(
    value = as.numeric(mean(r2)),
    elpd = as.numeric(loo$elpd %||% NA_real_)
  )
}

# Test-set metrics -------------------------------------------------------

#' @title ELPD Test-Set Metric
#' @description Expected log-predictive density on a held-out test set,
#'   estimated by log-sum-exp over posterior draws of `context$log_lik`.
#'   Returns NA when no test set.
#' @param name Character string naming the metric. Defaults to "elpd_test".
#' @return An `ElpdTestMetric` object.
#' @keywords internal
ElpdTestMetric <- S7::new_class(
  "ElpdTestMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "elpd_test"),
    needs = S7::new_property(S7::class_character, default = "log_lik"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname ElpdTestMetric
#' @description Constructor for ElpdTestMetric.
#' @return An `ElpdTestMetric` object.
#' @export
#' @examples
#' elpd_test_metric()
elpd_test_metric <- function(name = "elpd_test") {
  ElpdTestMetric(name = name, needs = "log_lik", required = FALSE)
}

S7::method(compute_metric, ElpdTestMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  test_data <- data_bundle$test
  if (is.null(test_data) || is.null(context$log_lik)) {
    return(list(value = NA_real_, n_obs = NA_integer_))
  }
  ll <- context$log_lik
  if (!is.matrix(ll)) return(list(value = NA_real_, n_obs = NA_integer_))
  elpd <- sum(apply(ll, 2, function(col) {
    m <- max(col)
    m + log(mean(exp(col - m)))
  }))
  list(value = elpd, n_obs = ncol(ll))
}

#' @title RMSE Test-Set Metric
#' @description Root-mean-square error of posterior-mean predictions on the
#'   held-out test set. Returns NA when no test set.
#' @param name Character string naming the metric. Defaults to "rmse_test".
#' @return A `RmseTestMetric` object.
#' @keywords internal
RmseTestMetric <- S7::new_class(
  "RmseTestMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "rmse_test"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname RmseTestMetric
#' @description Constructor for RmseTestMetric.
#' @return A `RmseTestMetric` object.
#' @export
#' @examples
#' rmse_test_metric()
rmse_test_metric <- function(name = "rmse_test") {
  RmseTestMetric(name = name, needs = "predictions", required = FALSE)
}

S7::method(compute_metric, RmseTestMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  test_data <- data_bundle$test
  if (is.null(test_data) || is.null(context$predictions)) {
    return(list(value = NA_real_, n_obs = NA_integer_))
  }
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean
  n <- min(length(actual), length(predicted))
  list(value = sqrt(mean((predicted[seq_len(n)] - actual[seq_len(n)])^2)), n_obs = n)
}

#' @title R-squared Test-Set Metric
#' @description Coefficient of determination of posterior-mean predictions on
#'   the held-out test set. Returns NA when no test set.
#' @param name Character string naming the metric. Defaults to "r2_test".
#' @return A `R2TestMetric` object.
#' @keywords internal
R2TestMetric <- S7::new_class(
  "R2TestMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "r2_test"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

#' @rdname R2TestMetric
#' @description Constructor for R2TestMetric.
#' @return A `R2TestMetric` object.
#' @export
#' @examples
#' r2_test_metric()
r2_test_metric <- function(name = "r2_test") {
  R2TestMetric(name = name, needs = "predictions", required = FALSE)
}

S7::method(compute_metric, R2TestMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  test_data <- data_bundle$test
  if (is.null(test_data) || is.null(context$predictions)) {
    return(list(value = NA_real_, n_obs = NA_integer_))
  }
  actual <- test_data[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean
  n <- min(length(actual), length(predicted))
  actual <- actual[seq_len(n)]
  predicted <- predicted[seq_len(n)]
  ss_res <- sum((actual - predicted)^2)
  ss_tot <- sum((actual - mean(actual))^2)
  list(value = 1 - ss_res / ss_tot, n_obs = n)
}
