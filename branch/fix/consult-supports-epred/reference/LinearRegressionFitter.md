# Conjugate Bayesian Linear Regression Fitter

Exact conjugate Normal-Inverse-Gamma (NIG) Bayesian linear regression.
Fits `y ~ N(X beta, sigma^2)` analytically and draws i.i.d. samples from
the joint posterior `(beta, sigma)`. No Stan, milliseconds per fit,
**real posteriors** — the package's teaching backbone (D1).

The model formula is taken from `fit_spec$formula` (a base R formula,
default `response ~ .`). Posterior draws use plain parameter names
(`Intercept`, `<coef>`, `sigma`) so they line up with
[`resolve_draw_columns()`](https://sims1253.github.io/bayesim/reference/resolve_draw_columns.md)
and the generators' cleaned names out of the box.

## Usage

``` r
LinearRegressionFitter(
  name = "linear_regression",
  supports_predictions = TRUE,
  supports_log_lik = TRUE,
  supports_loo = TRUE,
  supports_epred = TRUE,
  n_draws = 1000L,
  prior_mean = 0,
  prior_precision = 1e-06,
  a0 = 2,
  b0 = 1e-06
)
```

## Arguments

- name:

  Character string identifying the fitter.

- supports_predictions:

  Logical; whether predictions are supported.

- supports_log_lik:

  Logical; whether log-likelihood is supported.

- supports_loo:

  Logical; whether LOO-CV is supported.

- supports_epred:

  Logical; whether posterior expectation predictions
  ([`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md))
  are supported.

- n_draws:

  Positive integer; number of i.i.d. posterior draws.

- prior_mean:

  Numeric vector (length = number of coefficients, including intercept)
  or scalar; prior mean of `beta`. Recycled. Default 0.

- prior_precision:

  Numeric scalar; prior precision of `beta` per unit of `sigma^2` (i.e.
  `Lambda0 = prior_precision * I`). A small value gives a weak prior.
  Default `1e-6`.

- a0:

  Positive numeric; Inv-Gamma prior shape for `sigma^2`. Default `2`.

- b0:

  Positive numeric; Inv-Gamma prior rate for `sigma^2`. Default `1e-6`.

## Value

An S7 `LinearRegressionFitter` object.

## See also

[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[BrmsFitter](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fitter <- LinearRegressionFitter(n_draws = 500L)
} # }
```
