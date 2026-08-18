# Custom Fitters

## Custom Fitters

bayesim ships three fitters:
[`LinearRegressionFitter()`](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
(exact conjugate inference, no Stan),
[`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
(Stan via brms, with the compile-once model bank), and
[`CmdStanFitter()`](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md)
(your own `.stan` programs via cmdstanr). To use a different backend,
implement a custom fitter by extending the `Fitter` class. The core
contract is deliberately small: implement
[`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
and
[`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md),
then opt into optional capabilities by declaring the matching
`supports_*` property and method.

### The Fitter contract

The public fitter contract is:

| generic | signature | returns |
|----|----|----|
| `fit_model(fitter, data_bundle, fit_spec, seed, task_ctx)` | required | a `bayesim_fit_result`, constructed with exported [`new_fit_result()`](https://sims1253.github.io/bayesim/reference/new_fit_result.md) |
| `extract_draws(fitter, fit_result, variables)` | required | numeric draws matrix (rows = draws, cols = parameters, with unique colnames) |
| `predict_fit(fitter, fit_result, newdata, seed)` | optional; `supports_predictions = TRUE` | list with `predicted_mean` (length N), `predicted_samples` (S x N), `predicted_sd` (length N) |
| `log_lik_matrix(fitter, fit_result, newdata)` | optional; `supports_log_lik = TRUE` | numeric log-likelihood matrix (S x N) |
| `predict_epred(fitter, fit_result, newdata)` | optional; `supports_epred = TRUE` | numeric expectation matrix (S x N), or `NULL` when unsupported (`r2_loo` then degrades to NA) |
| `loo_fit(fitter, fit_result)` | optional; `supports_loo = TRUE` | list with `elpd`, `p_loo`, `elpd_se`, `pareto_k` |
| `fit_diagnostics(fitter, fit_result)` | optional | named diagnostic list; the default is [`list()`](https://rdrr.io/r/base/list.html) |

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

# fit_model(): estimate the model and package draws
S7::method(fit_model, LinearFitter) <- function(fitter, data_bundle, fit_spec, seed, task_ctx) {
  train <- data_bundle$train
  fit <- lm(y ~ x, data = train)
  coefs <- coef(fit)
  vc <- vcov(fit)
  draws <- MASS::mvrnorm(fitter@n_draws, mu = coefs, Sigma = vc)
  colnames(draws) <- c("intercept", "slope")

  # Package into the supported public fit-result type. The engine validates
  # and canonicalizes the draws at the task seam.
  bayesim::new_fit_result(
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
}

# extract_draws(): return the draws matrix
S7::method(extract_draws, LinearFitter) <- function(fitter, fit_result, variables = NULL) {
  fit_result$draws
}

# predict_fit(): posterior-mean predictions
# Convention: predicted_samples is S x N (draws as rows, observations as cols),
# matching log_lik() and predict_epred().
S7::method(predict_fit, LinearFitter) <- function(fitter, fit_result, newdata = NULL, seed = NULL) {
  data <- if (is.null(newdata)) fit_result$data_bundle$train else newdata
  # draws is S x P; design matrix X is N x P, so draws %*% t(X) is S x N.
  X <- as.matrix(cbind(1, data$x))
  preds <- fit_result$draws %*% t(X)
  list(
    predicted_mean = colMeans(preds),
    predicted_samples = preds,
    predicted_sd = apply(preds, 2, sd)
  )
}

# log_lik_matrix(): per-observation log-likelihood of each draw
S7::method(log_lik_matrix, LinearFitter) <- function(fitter, fit_result, newdata = NULL) {
  data <- if (is.null(newdata)) fit_result$data_bundle$train else newdata
  fit <- fit_result$fit
  sigma2 <- sum(residuals(fit)^2) / fit$df.residual
  X <- as.matrix(cbind(1, data$x))
  mu <- fit_result$draws %*% t(X)
  residual <- sweep(mu, 2, data$y, "-")
  # matrix: draws (rows) x observations (cols)
  -0.5 * log(2 * pi * sigma2) - residual^2 / (2 * sigma2)
}

# fit_diagnostics()
S7::method(fit_diagnostics, LinearFitter) <- function(fitter, fit_result) {
  fit_result$diagnostics
}
```

### Using the custom fitter

``` r

my_generator <- function(data_spec, task_ctx) {
  train_x <- stats::rnorm(data_spec$n_train)
  test_x <- stats::rnorm(data_spec$n_test)
  list(
    train = data.frame(
      y = data_spec$slope * train_x + stats::rnorm(data_spec$n_train),
      x = train_x
    ),
    test = data.frame(
      y = data_spec$slope * test_x + stats::rnorm(data_spec$n_test),
      x = test_x
    ),
    response = "y",
    true_params = c(slope = data_spec$slope),
    vars_of_interest = c("slope")
  )
}

# This metric deliberately requests both optional capabilities. Its output
# makes the S x N contract observable in the completed simulation.
ContractMetric <- S7::new_class(
  "ContractMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "contract"),
    needs = S7::new_property(
      S7::class_character,
      default = c("predictions", "log_lik")
    ),
    required = S7::new_property(S7::class_logical, default = TRUE)
  )
)

S7::method(compute_metric, ContractMetric) <- function(
  metric, fit_result, data_bundle, context, task_ctx
) {
  list(
    prediction_mean = mean(context$predictions$predicted_mean),
    mean_log_lik = mean(context$log_lik),
    prediction_n = ncol(context$predictions$predicted_samples),
    log_lik_n = ncol(context$log_lik)
  )
}

fitter <- LinearFitter()
representative_bundle <- withr::with_seed(
  11L,
  my_generator(list(n_train = 30L, n_test = 8L, slope = 2), list())
)
validate_fitter(
  fitter,
  smoke_test = TRUE,
  data_bundle = representative_bundle,
  fit_spec = list(model = "linear")
)

config <- simulation_config(
  data_grid = data.frame(n_train = 50L, n_test = 8L, slope = 2),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_generator,
  fitter = fitter,
  metrics = list(
    # Metric is abstract, so constructing ContractMetric() directly honors
    # the defaults declared above (name, needs = c("predictions",
    # "log_lik"), required = TRUE).
    ContractMetric(),
    posterior_summary_metric()
  ),
  n_replicates = 4L,
  seed = 42L
)

result <- run_simulation(config, progress = FALSE, verbose = FALSE)
stopifnot(
  all(result$summary$status == "success"),
  all(result$summary$contract__prediction_n == 8L),
  all(result$summary$contract__log_lik_n == 8L),
  all(is.finite(result$summary$contract__prediction_mean)),
  all(is.finite(result$summary$contract__mean_log_lik))
)
head(result$summary)
#>            task_id  status stop_reason contract__prediction_mean
#> 1 d001_f001_r00001 success        <NA>                -0.2697102
#> 2 d001_f001_r00002 success        <NA>                 0.5809961
#> 3 d001_f001_r00003 success        <NA>                -1.4781161
#> 4 d001_f001_r00004 success        <NA>                -0.7460078
#>   contract__mean_log_lik contract__prediction_n contract__log_lik_n
#> 1              -1.202266                      8                   8
#> 2              -1.255003                      8                   8
#> 3              -1.278065                      8                   8
#> 4              -1.591544                      8                   8
#>   posterior_summary__mean__slope posterior_summary__median__slope
#> 1                       1.828058                         1.820413
#> 2                       1.898103                         1.905611
#> 3                       1.990360                         1.997124
#> 4                       2.069987                         2.069358
#>   posterior_summary__sd__slope posterior_summary__q_lower__slope
#> 1                    0.1283973                          1.576287
#> 2                    0.1345020                          1.645288
#> 3                    0.1380528                          1.723293
#> 4                    0.1502059                          1.754710
#>   posterior_summary__q_upper__slope truth__slope rhat_max ess_bulk_min
#> 1                          2.071709            2        1          500
#> 2                          2.149415            2        1          500
#> 3                          2.259547            2        1          500
#> 4                          2.337851            2        1          500
#>   ess_tail_min divergent timing_total rep_idx data_n_train data_n_test
#> 1          500         0  0.014540195       1           50           8
#> 2          500         0  0.017550707       2           50           8
#> 3          500         0  0.002790213       3           50           8
#> 4          500         0  0.002588987       4           50           8
#>   data_slope fit_model
#> 1          2    linear
#> 2          2    linear
#> 3          2    linear
#> 4          2    linear
```

#### Parallel runs

`run_simulation(config, workers = N)` sets up mirai daemons for the run
and tears them down afterwards. `workers = 1` is *genuinely sequential*:
no daemons are launched and tasks run in-process. That matters for
fitters and metrics defined in *your* package or script —
[`S7::method()`](https://rconsortium.github.io/S7/reference/method.html)
registrations live in the process that ran them and do not travel to
daemon workers, so an external S7 fitter that passes
`validate_fitter(smoke_test = TRUE)` will lose
[`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
dispatch on daemons. Use `workers = 1` (or omit `workers`) for
externally defined fitters; parallel execution starts at `workers >= 2`.

### Validation

The executable workflow above calls `validate_fitter(smoke_test = TRUE)`
on representative data before running the study. That conformance check
executes fit, draw extraction, prediction, log-likelihood, and
diagnostics methods and checks their shapes. For a quick
declaration-only check, use:

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

See
[`?BrmsFitter`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
for the full method reference, and the `BrmsFitter` implementation in
`R/brms-fitter.R` of the bayesim source tree for a complete worked
example. Because the model bank keys on `fit_spec` columns (`formula`,
`family`, `prior`, `stanvars`), supplying these as `fit_grid` columns
lets the engine compile once and reuse the binary across all replicates
of the same condition.
