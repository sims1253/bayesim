# Mock Fitter for Testing (Internal)

A simple implementation of the Fitter class for **testing purposes
only**. This fitter returns predetermined values and does not perform
actual Bayesian inference.

**WARNING**: MockFitter ignores the `fit_spec` model specification. Do
NOT use for scientific simulation studies comparing different models.
For real inference, use
[BrmsFitter](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
or create a custom fitter.

For a lightweight no-Stan alternative, see the `LinearFitter` example in
[`vignette("custom-fitters")`](https://sims1253.github.io/bayesim/articles/custom-fitters.md).

## Usage

``` r
MockFitter(
  name = "mock",
  supports_predictions = TRUE,
  supports_log_lik = TRUE,
  supports_loo = TRUE,
  supports_epred = FALSE,
  n_draws = 100L,
  n_chains = 4L
)
```

## Arguments

- name:

  Character string identifying the fitter

- supports_predictions:

  Logical; whether predictions are supported

- supports_log_lik:

  Logical; whether log-likelihood is supported

- supports_loo:

  Logical; whether LOO-CV is supported

- n_draws:

  Integer; number of posterior draws to simulate

- n_chains:

  Integer; number of chains to simulate

## Value

An S7 class object representing a MockFitter

## See also

[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md) for the
abstract base class,
[BrmsFitter](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
and
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
for real inference
