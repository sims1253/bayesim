# Construct a fitter-agnostic prior-draws data generator

`prior_draws_generator()` is the fitter-agnostic analogue of
[`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md).
It works through the S7
[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md)
interface rather than brms-specific functions, so it can be used with
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md),
[BrmsFitter](https://sims1253.github.io/bayesim/reference/BrmsFitter.md),
[CmdStanFitter](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md),
or any custom Fitter.

## Usage

``` r
prior_draws_generator(
  fitter,
  fit_spec,
  pilot_bundle,
  predictor_generator,
  response = NULL,
  n_draws = NULL
)
```

## Arguments

- fitter:

  An S7 [Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md)
  object (e.g.
  [LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)).

- fit_spec:

  A list (single-row fit_grid entry) carrying at least `formula` (a base
  R formula). For
  [LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
  a formula like `y ~ x` is expected.

- pilot_bundle:

  A `data_bundle` list (`train`, `response`, etc.) used for the one-time
  preconditioning fit. The caller is responsible for providing a
  representative pilot dataset. Must contain a `train` data.frame whose
  column names match the design implied by `fit_spec$formula`.

- predictor_generator:

  Function `(data_spec, task_ctx) -> data.frame` producing the design
  matrix of predictors (everything except the response). Must consume
  the ambient RNG state.

- response:

  Name of the response column. Defaults to `pilot_bundle$response`,
  falling back to the LHS of `fit_spec$formula`, then to `"y"`.

- n_draws:

  Optional integer override for the number of draws to store; if `NULL`
  (default), uses the number of draws returned by
  [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md).

## Value

A generator function `(data_spec, task_ctx) -> data_bundle`.

## Details

The factory fits the model ONCE on `pilot_bundle` (provided by the
caller), extracts parameter draws via
[`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md),
and stores them. The returned closure, on each call, picks a draw
deterministically indexed by `task_ctx$rep_idx` (wrapped modulo the
number of stored draws), uses it as `true_params`, and forward-simulates
`y` from the supplied predictors.

Forward simulation: the response is drawn from a Gaussian with mean
equal to the linear predictor `X theta` (using the coefficient columns
of the draw) and standard deviation equal to the `sigma` column of the
draw (if present, else 1). This is the natural data-generating process
for Gaussian linear models — the common case for
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md).

## Limitations (non-brms)

A true model prior is directly accessible only for brms fits
([`brms::prior_draws()`](https://paulbuerkner.com/brms/reference/prior_draws.brmsfit.html)).
For other fitters there is no generic "prior_draws" S7 method, so this
factory falls back to the draws stored on the pilot fit. For
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
those are NIG-prior-conditioned posterior draws (the prior is weak by
default, `prior_precision = 1e-6`), which makes this an *approximate*
prior-predictive path concentrated on the pilot's posterior region. Brms
users who need full prior-predictive coverage should prefer the
brms-specific
[`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md).

## See also

[`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md),
[`forward_sim_generator()`](https://sims1253.github.io/bayesim/reference/forward_sim_generator.md),
[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
