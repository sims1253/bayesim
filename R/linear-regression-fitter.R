# LinearRegressionFitter --------------------------------------------------
#
# An exact conjugate Normal-Inverse-Gamma (NIG) Bayesian linear regression
# fitter. This is bayesim's executable-docs teaching backbone: students get
# real posteriors, real coverage ~0.95, real SBC uniformity, in milliseconds
# and with zero Stan. See D1.
#
# Model: y ~ N(X beta, sigma^2), with NIG prior
#   beta | sigma^2 ~ N(mu0, sigma^2 * Lambda0^{-1})
#   sigma^2        ~ Inv-Gamma(a0, b0)
# Posterior (conjugate):
#   Lambda_n = Lambda0 + X'X
#   mu_n     = Lambda_n^{-1} (Lambda0 mu0 + X' y)
#   a_n      = a0 + n/2
#   b_n      = b0 + 0.5 (y' y + mu0' Lambda0 mu0 - mu_n' Lambda_n mu_n)
# Draws: sigma^2 ~ Inv-Gamma(a_n, b_n); beta | sigma^2 ~ N(mu_n, sigma^2 Lambda_n^{-1}).
# All matrices follow the draws-x-observations (S x N) convention.

#' @title Conjugate Bayesian Linear Regression Fitter
#' @description
#' Exact conjugate Normal-Inverse-Gamma (NIG) Bayesian linear regression. Fits
#' `y ~ N(X beta, sigma^2)` analytically and draws i.i.d. samples from the
#' joint posterior `(beta, sigma)`. No Stan, milliseconds per fit, **real
#' posteriors** — the package's teaching backbone (D1).
#'
#' The model formula is taken from `fit_spec$formula` (a base R formula,
#' default `response ~ .`). Posterior draws use plain parameter names
#' (`Intercept`, `<coef>`, `sigma`) so they line up with
#' [resolve_draw_columns()] and the generators' cleaned names out of the box.
#'
#' @param name Character string identifying the fitter.
#' @param supports_predictions Logical; whether predictions are supported.
#' @param supports_log_lik Logical; whether log-likelihood is supported.
#' @param supports_loo Logical; whether LOO-CV is supported.
#' @param supports_epred Logical; whether posterior expectation predictions
#'   (`predict_epred()`) are supported.
#' @param n_draws Positive integer; number of i.i.d. posterior draws.
#' @param prior_mean Numeric vector (length = number of coefficients, including
#'   intercept) or scalar; prior mean of `beta`. Recycled. Default 0.
#' @param prior_precision Numeric scalar; prior precision of `beta` per unit of
#'   `sigma^2` (i.e. `Lambda0 = prior_precision * I`). A small value gives a
#'   weak prior. Default `1e-6`.
#' @param a0 Positive numeric; Inv-Gamma prior shape for `sigma^2`. Default `2`.
#' @param b0 Positive numeric; Inv-Gamma prior rate for `sigma^2`. Default `1e-6`.
#'
#' @return An S7 `LinearRegressionFitter` object.
#' @export
#' @seealso [Fitter], [BrmsFitter]
#' @examples
#' \dontrun{
#' fitter <- LinearRegressionFitter(n_draws = 500L)
#' }
LinearRegressionFitter <- S7::new_class(
  "LinearRegressionFitter",
  parent = Fitter,
  properties = list(
    name = S7::new_property(S7::class_character, default = "linear_regression"),
    supports_predictions = S7::new_property(S7::class_logical, default = TRUE),
    supports_log_lik = S7::new_property(S7::class_logical, default = TRUE),
    supports_loo = S7::new_property(S7::class_logical, default = TRUE),
    supports_epred = S7::new_property(S7::class_logical, default = TRUE),
    n_draws = S7::new_property(S7::class_integer, default = 1000L),
    prior_mean = S7::new_property(
      S7::new_union(S7::class_numeric, NULL),
      default = 0
    ),
    prior_precision = S7::new_property(S7::class_numeric, default = 1e-6),
    a0 = S7::new_property(S7::class_numeric, default = 2),
    b0 = S7::new_property(S7::class_numeric, default = 1e-6)
  )
)

# Internal: build the design matrix and response vector from a data.frame and a
# formula. Returns list(X, y, coef_names). The intercept column is named
# "Intercept" (plain, matching brms' cleaned names).
.lm_design <- function(data, response, formula) {
  if (!is.null(formula)) {
    mf <- stats::model.frame(formula, data = data, na.action = na.pass)
    X <- stats::model.matrix(formula, mf)
    y <- stats::model.response(mf)
  } else {
    # Default response ~ .
    predictors <- setdiff(names(data), response)
    if (length(predictors) == 0L) {
      stop(bayesim_contract_error(
        "LinearRegressionFitter requires at least one predictor column"
      ))
    }
    y <- data[[response]]
    # model.matrix with intercept on the predictor columns
    fml <- stats::as.formula(paste("~", paste(predictors, collapse = " + ")))
    X <- stats::model.matrix(fml, data)
  }
  # Plain coefficient names: "Intercept" instead of "(Intercept)".
  cn <- colnames(X)
  cn[cn == "(Intercept)"] <- "Intercept"
  colnames(X) <- cn

  # Missing values would silently propagate through model.matrix into NaN
  # posteriors; reject them up front as a recoverable data error.
  if (anyNA(X) || anyNA(y)) {
    stop(bayesim_data_error(paste(
      "LinearRegressionFitter received data with missing values;",
      "NA predictors or responses are not supported (posterior would be NaN)."
    )))
  }

  list(X = X, y = as.numeric(y), coef_names = cn)
}

# Internal: compute the NIG posterior parameters.
.nig_posterior <- function(X, y, mu0, lambda0_scalar, a0, b0) {
  n <- nrow(X)
  p <- ncol(X)
  Lambda0 <- diag(lambda0_scalar, p)
  mu0 <- rep(mu0, length.out = p)
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  Lambda_n <- Lambda0 + XtX
  Lambda_n_inv <- solve(Lambda_n)
  mu_n <- drop(Lambda_n_inv %*% (Lambda0 %*% mu0 + Xty))
  a_n <- a0 + n / 2
  b_n <- b0 +
    0.5 *
      (sum(y^2) +
        crossprod(mu0, Lambda0 %*% mu0) -
        crossprod(mu_n, Lambda_n %*% mu_n))
  list(
    mu_n = mu_n,
    Lambda_n_inv = Lambda_n_inv,
    a_n = a_n,
    b_n = b_n,
    n = n,
    p = p
  )
}

# fit_model ----------------------------------------------------------------
S7::method(fit_model, LinearRegressionFitter) <- function(
  fitter,
  data_bundle,
  fit_spec,
  seed,
  task_ctx
) {
  timer <- make_timer()
  timer$start()

  data <- data_bundle$train
  if (is.null(data)) {
    stop(bayesim_contract_error(
      "LinearRegressionFitter requires data_bundle$train"
    ))
  }
  response <- data_bundle$response %||%
    all.vars(fit_spec$formula)[1L] %||%
    "y"
  formula <- fit_spec$formula

  des <- .lm_design(data, response, formula)
  post <- .nig_posterior(
    des$X,
    des$y,
    fitter@prior_mean,
    fitter@prior_precision,
    fitter@a0,
    fitter@b0
  )

  n_draws <- as.integer(fitter@n_draws)
  draws <- withr::with_seed(
    seed %||% 0L,
    .nig_draws(post, n_draws, des$coef_names)
  )

  timer$stop()
  # Store everything needed by the other methods in `fit`.
  fit_obj <- list(
    fitter = "linear_regression",
    data_bundle = data_bundle,
    design = des,
    posterior = post,
    formula = formula,
    response = response,
    seed = seed %||% 0L,
    n_draws = n_draws
  )

  new_fit_result(
    success = TRUE,
    fit = fit_obj,
    draws = draws,
    diagnostics = list(
      rhat_max = 1,
      ess_bulk = n_draws,
      ess_tail = n_draws,
      divergent = 0L,
      max_treedepth = 0L
    ),
    timing = list(
      total = timer$elapsed(),
      warmup = 0,
      sample = timer$elapsed()
    ),
    warnings = character(),
    error = NULL,
    data_bundle = data_bundle
  )
}

# Internal: i.i.d. joint posterior draws. Returns S x P matrix with colnames.
.nig_draws <- function(post, n_draws, coef_names) {
  # sigma^2 ~ Inv-Gamma(a_n, b_n): draw via 1 / Gamma(a_n, rate=b_n)
  sigma2 <- 1 / stats::rgamma(n_draws, shape = post$a_n, rate = post$b_n)
  # beta | sigma^2 ~ N(mu_n, sigma^2 * Lambda_n^{-1})
  p <- post$p
  beta_draws <- matrix(NA_real_, nrow = n_draws, ncol = p)
  chol_inv <- t(chol(post$Lambda_n_inv))
  for (s in seq_len(n_draws)) {
    beta_draws[s, ] <- post$mu_n +
      sqrt(sigma2[s]) * (chol_inv %*% stats::rnorm(p))
  }
  colnames(beta_draws) <- coef_names
  # Append sigma column.
  out <- cbind(beta_draws, sigma = sqrt(sigma2))
  out
}

# extract_draws -----------------------------------------------------------
S7::method(extract_draws, LinearRegressionFitter) <- function(
  fitter,
  fit_result,
  variables = NULL
) {
  draws <- if (inherits(fit_result, "bayesim_fit_result")) {
    fit_result$draws
  } else {
    fit_result
  }
  if (is.null(draws)) {
    return(NULL)
  }
  if (!is.null(variables)) {
    draws <- draws[, variables, drop = FALSE]
  }
  draws
}

# predict_epred: X %*% t(beta_draws) => S x N --------------------------------
S7::method(predict_epred, LinearRegressionFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL
) {
  fit_obj <- fit_result$fit
  data <- newdata %||%
    fit_result$data_bundle$train %||%
    fit_obj$data_bundle$train
  des <- .lm_design(data, fit_obj$response, fit_obj$formula)
  beta_draws <- fit_result$draws[, fit_obj$design$coef_names, drop = FALSE]
  # S x N: (S x P) %*% t(N x P)
  beta_draws %*% t(des$X)
}

# predict_fit: epred + noise -----------------------------------------------
S7::method(predict_fit, LinearRegressionFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL,
  seed = NULL
) {
  fit_obj <- fit_result$fit
  data <- newdata %||%
    fit_result$data_bundle$train %||%
    fit_obj$data_bundle$train
  des <- .lm_design(data, fit_obj$response, fit_obj$formula)
  draws <- fit_result$draws
  beta_draws <- draws[, fit_obj$design$coef_names, drop = FALSE]
  sigma_draws <- draws[, "sigma"]
  n_obs <- nrow(des$X)
  draw_predictions <- function() {
    mu <- beta_draws %*% t(des$X) # S x N
    noise <- matrix(
      stats::rnorm(n_obs * nrow(mu)),
      nrow = nrow(mu),
      ncol = n_obs
    )
    predicted_samples <- mu + sigma_draws * noise
    list(
      predicted_mean = drop(colMeans(predicted_samples)),
      predicted_samples = predicted_samples,
      predicted_sd = apply(predicted_samples, 2, stats::sd)
    )
  }
  # seed = NULL consumes the ambient RNG stream (required by the
  # forward-simulation generators, which need fresh noise per task); an
  # explicit seed gives reproducible predictions without touching that stream.
  if (is.null(seed)) {
    draw_predictions()
  } else {
    withr::with_seed(seed, draw_predictions())
  }
}

# log_lik_matrix: exact normal density, S x N ------------------------------
S7::method(log_lik_matrix, LinearRegressionFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL
) {
  fit_obj <- fit_result$fit
  data <- newdata %||%
    fit_result$data_bundle$train %||%
    fit_obj$data_bundle$train
  des <- .lm_design(data, fit_obj$response, fit_obj$formula)
  draws <- fit_result$draws
  beta_draws <- draws[, fit_obj$design$coef_names, drop = FALSE]
  sigma_draws <- draws[, "sigma"]
  mu <- beta_draws %*% t(des$X) # S x N
  resid <- sweep(mu, 2, des$y, "-") # subtract y from each column
  # log N(y | mu, sigma^2); sigma_draws is length S, broadcast over columns
  ll <- -0.5 * log(2 * pi) - log(sigma_draws) - 0.5 * (resid / sigma_draws)^2
  ll
}

# loo_fit ------------------------------------------------------------------
S7::method(loo_fit, LinearRegressionFitter) <- function(
  fitter,
  fit_result,
  log_lik = NULL
) {
  # #73: reuse a caller-supplied train-set matrix; standalone callers pass
  # NULL and pay the single computation here.
  ll <- log_lik %||% log_lik_matrix(fitter, fit_result)
  if (is.null(ll)) {
    return(list(
      elpd = NA_real_,
      p_loo = NA_real_,
      elpd_se = NA_real_,
      pareto_k = numeric(),
      r_eff = NULL
    ))
  }
  # i.i.d. draws => no chains; relative_eff is degenerate. PSIS-LOO still valid.
  loo_result <- tryCatch(
    suppressWarnings(loo::loo(ll)),
    error = function(e) NULL
  )
  if (is.null(loo_result)) {
    return(list(
      elpd = NA_real_,
      p_loo = NA_real_,
      elpd_se = NA_real_,
      pareto_k = numeric(),
      r_eff = NULL
    ))
  }
  list(
    elpd = loo_result$estimates["elpd_loo", "Estimate"],
    p_loo = loo_result$estimates["p_loo", "Estimate"],
    elpd_se = loo_result$estimates["elpd_loo", "SE"],
    pareto_k = loo::pareto_k_values(loo_result),
    r_eff = NULL
  )
}

# fit_diagnostics ----------------------------------------------------------
S7::method(fit_diagnostics, LinearRegressionFitter) <- function(
  fitter,
  fit_result
) {
  n_draws <- if (!is.null(fit_result$draws)) {
    nrow(fit_result$draws)
  } else {
    NA_integer_
  }
  list(
    rhat_max = 1,
    ess_bulk = n_draws,
    ess_tail = n_draws,
    divergent = 0L,
    max_treedepth = 0L
  )
}
