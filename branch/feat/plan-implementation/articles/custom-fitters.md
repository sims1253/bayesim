# Custom Fitters

## Custom Fitters

bayesim ships three fitters:
[`LinearRegressionFitter()`](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
(exact conjugate inference, no Stan),
[`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
(Stan via brms, with the compile-once model bank), and
[`CmdStanFitter()`](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md)
(your own `.stan` programs via cmdstanr). To use a different backend,
implement a custom fitter by extending the `Fitter` class and providing
methods for the required generics.

### The Fitter contract

Every fitter must implement these S7 generics:

| generic | signature | returns |
|----|----|----|
| `fit_model(fitter, data_bundle, fit_spec, seed, task_ctx)` | training data + spec | a `bayesim_fit_result` (use [`new_fit_result()`](https://sims1253.github.io/bayesim/reference/new_fit_result.md) internally) |
| `extract_draws(fitter, fit_result, variables)` | a fit result | a draws matrix (rows = draws, cols = parameters, with colnames) |
| `predict_fit(fitter, fit_result, newdata, seed)` | a fit result | list with `predicted_mean` (length N), `predicted_samples` (S x N), `predicted_sd` (length N) |
| `log_lik_matrix(fitter, fit_result, newdata)` | a fit result | log-likelihood matrix (S x N) |
| `loo_fit(fitter, fit_result)` | a fit result | list with `elpd`, `p_loo`, `elpd_se`, `pareto_k` |
| `fit_diagnostics(fitter, fit_result)` | a fit result | list with `rhat_max`, `ess_bulk_min`, `divergent`, … |

**All matrices are draws x observations (S x N)** — `predicted_samples`,
[`log_lik_matrix()`](https://sims1253.github.io/bayesim/reference/log_lik_matrix.md),
and
[`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md)
all put posterior draws in rows and observations in columns, matching
the brms/loo convention. This is the single most common custom-fitter
bug; `validate_fitter(smoke_test = TRUE)` rejects transposed matrices.

The generics avoid masking common names: `loo_fit` (not
[`loo::loo`](https://mc-stan.org/loo/reference/loo.html)),
`log_lik_matrix` (not
[`brms::log_lik`](https://mc-stan.org/rstantools/reference/log_lik.html)),
`fit_model` (not
[`generics::fit`](https://generics.r-lib.org/reference/fit.html)).

### A minimal fitter: LinearFitter

This fitter fits a linear model by OLS and synthesizes posterior draws
by bootstrapping the residual variance. It is executable without Stan.

``` r

library(bayesim)

# Define the class
LinearFitter <- S7::new_class(
  "LinearFitter",
  parent = Fitter,
  properties = list(
    name = S7::new_property(S7::class_character, default = "linear"),
    supports_predictions = S7::new_property(S7::class_logical, default = TRUE),
    supports_log_lik = S7::new_property(S7::class_logical, default = TRUE),
    supports_loo = S7::new_property(S7::class_logical, default = FALSE),
    n_draws = S7::new_property(S7::class_integer, default = 500L)
  )
)

# fit(): estimate the model and package draws
S7::method(fit_model, LinearFitter) <- function(fitter, data_bundle, fit_spec, seed, task_ctx) {
  train <- data_bundle$train
  fit <- lm(y ~ x, data = train)
  coefs <- coef(fit)
  vc <- vcov(fit)
  draws <- MASS::mvrnorm(fitter@n_draws, mu = coefs, Sigma = vc)
  colnames(draws) <- c("intercept", "slope")

  # Package into a fit result. new_fit_result() is an internal helper.
  result <- list(
    fit = fit,
    draws = draws,
    diagnostics = list(
      rhat_max = 1.0, ess_bulk_min = as.numeric(fitter@n_draws),
      ess_tail_min = as.numeric(fitter@n_draws), divergent = 0L
    ),
    timing = list(total = 0, warmup = 0, sample = 0),
    success = TRUE,
    warnings = character(0),
    error = NULL,
    data_bundle = data_bundle
  )
  class(result) <- "bayesim_fit_result"
  result
}

# extract_draws(): return the draws matrix
S7::method(extract_draws, LinearFitter) <- function(fitter, fit_result, variables = NULL) {
  fit_result$draws
}

# predict_fit(): posterior-mean predictions
# Convention: predicted_samples is S x N (draws as rows, observations as cols),
# matching log_lik() and predict_epred().
S7::method(predict_fit, LinearFitter) <- function(fitter, fit_result, newdata = NULL, seed = NULL) {
  data <- fit_result$data_bundle$train
  # draws is S x P; design matrix X is N x P, so draws %*% t(X) is S x N.
  X <- as.matrix(cbind(1, data$x))
  preds <- fit_result$draws %*% t(X)
  list(
    predicted_mean = colMeans(preds),
    predicted_samples = preds,
    predicted_sd = apply(preds, 2, sd)
  )
}

# log_lik(): per-observation log-likelihood of each draw
S7::method(log_lik_matrix, LinearFitter) <- function(fitter, fit_result, newdata = NULL) {
  data <- newdata %||% fit_result$data_bundle$train
  fit <- fit_result$fit
  sigma2 <- sum(residuals(fit)^2) / fit$df.residual
  X <- as.matrix(cbind(1, data$x))
  mu <- t(fit_result$draws %*% t(X))
  # matrix: draws (rows) x observations (cols)
  t(dnorm(t(data$y), mean = mu, sd = sqrt(sigma2), log = TRUE))
}

# loo_fit(): not supported by this fitter (supports_loo = FALSE)
S7::method(loo_fit, LinearFitter) <- function(fitter, fit_result) {
  list(elpd = NA_real_, p_loo = NA_real_, elpd_se = NA_real_, pareto_k = NA_real_)
}

# diagnostics()
S7::method(fit_diagnostics, LinearFitter) <- function(fitter, fit_result) {
  fit_result$diagnostics
}
```

### Using the custom fitter

``` r

my_generator <- function(data_spec, task_ctx) {
  n <- data_spec$n
  x <- stats::rnorm(n)
  y <- data_spec$slope * x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x), test = NULL, response = "y",
    true_params = c(slope = data_spec$slope),
    vars_of_interest = c("slope")
  )
}

config <- simulation_config(
  data_grid = data.frame(n = 50, slope = 2),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_generator,
  fitter = LinearFitter(),
  metrics = list(pred_rmse_metric(), pred_bias_metric(), posterior_summary_metric()),
  n_replicates = 4L,
  seed = 42L
)

result <- run_simulation(config, progress = FALSE)
#> 4 tasks = 1 data x 1 fit x 4 reps
#> ℹ Starting simulation with 4 tasks
#> Warning: pred_rmse_metric: no test set in data_bundle; returning NA.
#> ℹ Provide a test set via the data generator to measure predictive error.
#> Warning: pred_bias_metric: no test set in data_bundle; returning NA.
#> ℹ Provide a test set via the data generator to measure predictive error.
head(result$summary)
#>            task_id  status rmse__value rmse__n_obs bias__value
#> 1 d001_f001_r00001 success          NA          NA          NA
#> 2 d001_f001_r00002 success          NA          NA          NA
#> 3 d001_f001_r00003 success          NA          NA          NA
#> 4 d001_f001_r00004 success          NA          NA          NA
#>   posterior_summary__mean__slope posterior_summary__median__slope
#> 1                       1.855743                         1.848816
#> 2                       1.851992                         1.860257
#> 3                       2.086556                         2.094179
#> 4                       1.766086                         1.766338
#>   posterior_summary__sd__slope posterior_summary__q_lower__slope
#> 1                    0.1313655                          1.597739
#> 2                    0.1334661                          1.599657
#> 3                    0.1513471                          1.791782
#> 4                    0.1466753                          1.458247
#>   posterior_summary__q_upper__slope truth__slope rhat_max ess_bulk_min
#> 1                          2.106035            2        1          500
#> 2                          2.096462            2        1          500
#> 3                          2.382663            2        1          500
#> 4                          2.027047            2        1          500
#>   ess_tail_min divergent timing_total rep_idx data_n data_slope fit_model
#> 1          500         0  0.068969011       1     50          2    linear
#> 2          500         0  0.024803162       2     50          2    linear
#> 3          500         0  0.003398418       3     50          2    linear
#> 4          500         0  0.003047466       4     50          2    linear
```

### Validation

Use
[`validate_fitter()`](https://sims1253.github.io/bayesim/reference/validate_fitter.md)
to confirm a fitter implements all required generics before running a
study:

``` r

validate_fitter(LinearFitter())
```

## Raw Stan: CmdStanFitter

When brms cannot express your model, write the `.stan` file yourself and
wire it in with
[`CmdStanFitter()`](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md)
— no custom fitter class needed. You supply the Stan program, a
`stan_data` function mapping a data bundle to the Stan `data` list, and
(optionally) the names of `log_lik` / `epred` generated quantities:

``` r

stan_file <- system.file("stan", "linear_regression.stan", package = "bayesim")

fitter <- CmdStanFitter(
  stan_file = stan_file,
  stan_data = function(data_bundle, fit_spec) {
    X <- stats::model.matrix(~x, data_bundle$train)
    list(
      N = nrow(X), K = ncol(X) - 1L, X = X,
      y = data_bundle$train[[data_bundle$response]]
    )
  },
  log_lik = "log_lik",   # generated-quantities vector -> LOO metrics work
  epred = "mu",          # optional: in-sample expectation matrix
  chains = 2L, iter_warmup = 500L, iter_sampling = 500L
)
```

The shipped `inst/stan/linear_regression.stan` shows the pattern: a
`generated quantities` block computing `log_lik[n]` (pointwise
log-likelihood) and `mu[n]` (the expectation). With `log_lik` declared,
`supports_log_lik`/`supports_loo` are set automatically and
[`elpd_loo_metric()`](https://sims1253.github.io/bayesim/reference/ElpdLooMetric.md)
works out of the box; without it, an informative error names the
convention.

Limitations, by design: raw Stan has no newdata semantics, so
`supports_predictions` is `FALSE` (test-set prediction metrics need brms
or a custom fitter), and
[`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md)
returns the in-sample GQ matrix only when `epred` is declared.
Compilation is cached by cmdstanr (keyed on the file hash), so daemons
compile-or-cache-hit on first use — no model-bank involvement.

## Custom Metrics

Metrics are the other half of the extension surface — see
[`vignette("custom-metrics")`](https://sims1253.github.io/bayesim/articles/custom-metrics.md)
for the Metric contract, the output schema, `summary_type`, and
externalization of large outputs.

## A brms-based fitter (sketch)

For a real Stan fitter,
[`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
is the standard path and benefits from the model bank (one compile per
distinct spec). A custom brms fitter would:

1.  Call
    [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
    (or reuse a prefit via the model bank) in
    [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md).
2.  Use
    [`brms::as_draws_matrix()`](https://mc-stan.org/posterior/reference/draws_matrix.html)
    in
    [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md).
3.  Use
    [`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
    /
    [`brms::log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
    in the other generics.

See the [`BrmsFitter`
source](https://github.com/sims1253/bayesim/blob/main/R/brms-fitter.R)
for a complete reference implementation. Because the model bank keys on
`fit_spec` columns (`formula`, `family`, `prior`, `stanvars`), supplying
these as `fit_grid` columns lets the engine compile once and reuse the
binary across all replicates of the same condition.
