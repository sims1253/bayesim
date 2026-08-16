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
  supports_predictions = FALSE,
  supports_log_lik = FALSE,
  supports_loo = FALSE,
  supports_epred = FALSE
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

- supports_epred:

  Logical indicating if the fitter supports posterior expectation
  predictions
  ([`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md);
  required by the `r2_loo` / `rmse_loo` LOO metrics)

## Value

An S7 class object representing the abstract Fitter

## Methods

The following S7 generics form the fitter interface. A minimal custom
fitter only needs to implement
[`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
and
[`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md);
diagnostics default to an empty list and unsupported optional
capabilities default to `NULL`.

- `fit_model(fitter, data_bundle, fit_spec, seed, task_ctx)`:

  Main fitting method

- `extract_draws(fitter, fit_result, variables = NULL)`:

  Extract posterior draws

- `predict_fit(fitter, fit_result, newdata = NULL, seed = NULL)`:

  Generate predictions

- `log_lik_matrix(fitter, fit_result, newdata = NULL)`:

  Pointwise log-likelihood

- `loo_fit(fitter, fit_result)`:

  LOO-CV computation

- `fit_diagnostics(fitter, fit_result)`:

  Extract fit diagnostics

## Creating Custom Fitters

To create a custom fitter, extend this class and implement methods for
the core S7 generics:
[`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
and
[`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md).
Implement optional
[`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md),
[`log_lik_matrix()`](https://sims1253.github.io/bayesim/reference/log_lik_matrix.md),
and
[`loo_fit()`](https://sims1253.github.io/bayesim/reference/loo_fit.md)
methods only when the matching `supports_*` property is `TRUE`. All
matrices follow the draws-by-observations (S x N) orientation; see
[`vignette("custom-fitters")`](https://sims1253.github.io/bayesim/articles/custom-fitters.md)
for the full contract.

## See also

[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md),
[BrmsFitter](https://sims1253.github.io/bayesim/reference/BrmsFitter.md),
and
[CmdStanFitter](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md)
for the built-in implementations, and
[`validate_fitter()`](https://sims1253.github.io/bayesim/reference/validate_fitter.md)
to check a custom fitter against the contract
