# Abstract Fitter Class

Abstract base class for Bayesian model fitters in bayesim.

The Fitter class defines the interface that all model fitters must
implement. It provides a consistent API for fitting Bayesian models,
extracting posterior draws, generating predictions, computing
log-likelihoods, and performing model diagnostics.

## Usage

``` r
Fitter(
  name = character(0),
  supports_predictions = TRUE,
  supports_log_lik = TRUE,
  supports_loo = TRUE
)
```

## Arguments

- name:

  Character string identifying the fitter (e.g., "stan", "brms")

- supports_predictions:

  Logical indicating if the fitter supports predictions

- supports_log_lik:

  Logical indicating if the fitter supports log-likelihood computation

- supports_loo:

  Logical indicating if the fitter supports LOO-CV

## Value

An S7 class object representing the abstract Fitter

## Methods

The following S7 generics must be implemented by subclasses:

- `fit(fitter, data_bundle, fit_spec, seed, task_ctx)`:

  Main fitting method

- `extract_draws(fitter, fit_result, variables = NULL)`:

  Extract posterior draws

- `predict_fit(fitter, fit_result, newdata = NULL, seed = NULL)`:

  Generate predictions

- `log_lik(fitter, fit_result, newdata = NULL)`:

  Pointwise log-likelihood

- `loo(fitter, fit_result)`:

  LOO-CV computation

- `diagnostics(fitter, fit_result)`:

  Extract fit diagnostics

## Creating Custom Fitters

To create a custom fitter, extend this class and implement methods for
the S7 generics:
[`fit()`](https://sims1253.github.io/bayesim/reference/fit.md),
[`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`log_lik()`](https://sims1253.github.io/bayesim/reference/log_lik.md),
[`loo()`](https://sims1253.github.io/bayesim/reference/loo.md),
[`diagnostics()`](https://sims1253.github.io/bayesim/reference/diagnostics.md).

## See also

[MockFitter](https://sims1253.github.io/bayesim/reference/MockFitter.md)
for a simple example implementation
