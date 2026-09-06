# Construct a prior-predictive data generator

Draws a parameter vector theta from the model prior (via a
`sample_prior = "only"` brmsfit) and simulates data `y ~ p(y | theta)`
using
[`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html).
Each task uses a distinct prior draw, deterministically indexed by
`task_ctx$rep_idx`, so Simulation-Based Calibration ranks are
well-defined and resume is reproducible.

## Usage

``` r
prior_predictive_generator(
  prior_fit,
  predictor_generator = NULL,
  vars_of_interest = NULL,
  response = NULL
)
```

## Arguments

- prior_fit:

  A brmsfit compiled with `sample_prior = "only"` (or a formula +
  family + prior combination to be compiled; see Details). Must contain
  prior draws.

- predictor_generator:

  Function `(data_spec, task_ctx) -> data.frame` producing the design
  matrix of predictors (everything except the response). Must consume
  the ambient RNG state. If `NULL`, the prior_fit's own data is reused
  at its original size.

- vars_of_interest:

  Character vector naming the prior parameters to report as
  `true_params` (defaults to all population-level effects `"b_<name>"`,
  renamed to `<name>`).

- response:

  Name of the response column (defaults to the LHS of `prior_fit`'s
  formula).

## Value

A generator function `(data_spec, task_ctx) -> data_bundle`.

## Details

The prior model is compiled once (a `sample_prior = "only"` brmsfit);
reuse it across tasks via the model bank or by passing the same object.
Predictor covariates not implied by the prior are supplied by
`predictor_generator`.
