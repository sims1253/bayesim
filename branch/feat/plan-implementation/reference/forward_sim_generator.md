# Construct a fitter-agnostic forward-simulation (IFS) data generator

`forward_sim_generator()` is the fitter-agnostic analogue of
[`ifs_generator()`](https://sims1253.github.io/bayesim/reference/ifs_generator.md).
It fits the model ONCE on `pilot_bundle`, draws theta from the posterior
(via
[`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)),
and forward-simulates `y` via the fitter's
[`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
(posterior-predictive) at the chosen draw's parameters. Works with any
Fitter that supports
[`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
(LinearRegressionFitter, BrmsFitter, ...).

## Usage

``` r
forward_sim_generator(
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
  object supporting
  [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md).

- fit_spec:

  A list (single-row fit_grid entry) with at least `formula`.

- pilot_bundle:

  A `data_bundle` list used for the one-time preconditioning fit.

- predictor_generator:

  Function `(data_spec, task_ctx) -> data.frame` producing predictor
  covariates. Must consume the ambient RNG state.

- response:

  Name of the response column. Defaults to `pilot_bundle$response`,
  falling back to the LHS of `fit_spec$formula`, then to `"y"`.

- n_draws:

  Optional integer override for the number of draws to store.

## Value

A generator function `(data_spec, task_ctx) -> data_bundle`.

## Details

Unlike
[`prior_draws_generator()`](https://sims1253.github.io/bayesim/reference/prior_draws_generator.md)
(which targets the prior), this generator concentrates the truth draw in
a region of high posterior mass — the canonical SBC generator for models
with diffuse or improper priors.

Because forward simulation here relies on
[`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md),
the response is drawn exactly as the fitter implements its
posterior-predictive sampling, which respects the fitter's response
distribution and link function. Each task uses a distinct draw,
deterministically indexed by `task_ctx$rep_idx`.

## Limitations

The `true_params` reported by this generator are the
[`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)
columns for the selected draw; for
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
these are `Intercept`, `<coef>`, `sigma`. The response is
forward-simulated via the fitter's
[`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
applied to the predictor design (a single Gaussian draw for Gaussian
fitters). Fitters without
[`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
support are not supported (e.g. raw
[CmdStanFitter](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md),
which has no newdata semantics).

## See also

[`ifs_generator()`](https://sims1253.github.io/bayesim/reference/ifs_generator.md),
[`prior_draws_generator()`](https://sims1253.github.io/bayesim/reference/prior_draws_generator.md),
[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
