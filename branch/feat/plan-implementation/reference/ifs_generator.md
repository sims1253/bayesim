# Construct an inverse forward sampling (IFS) data generator

Draws a parameter vector theta from a preconditioning fit (typically a
posterior fit on a pilot dataset), uses it as the data-generating truth,
and forward-simulates data `y ~ p(y | theta)` respecting multivariate
response dependencies. Each task uses a distinct draw, deterministically
indexed by `task_ctx$rep_idx`.

## Usage

``` r
ifs_generator(
  prefit,
  predictor_generator = NULL,
  vars_of_interest = NULL,
  response = NULL,
  lower_bound = NULL,
  upper_bound = NULL,
  truncate = FALSE
)
```

## Arguments

- prefit:

  A brmsfit with posterior draws to sample theta from (the
  preconditioning fit).

- predictor_generator:

  Function `(data_spec, task_ctx) -> data.frame` producing predictor
  covariates. Must consume the ambient RNG state. If `NULL`,
  `prefit$data` is reused.

- vars_of_interest:

  Character vector naming the parameters to report as `true_params`.
  Defaults to population-level effects.

- response:

  Name of the response column; defaults to the LHS of `prefit`'s
  formula.

- lower_bound, upper_bound:

  Optional numeric bounds for the response domain. If both `NULL`
  (default), no bounds are applied. If set, the out-of-bounds policy is
  governed by `truncate`:

  - `truncate = FALSE` (default): out-of-bounds draws set the response
    to all `NA`, which fails downstream validation and drops the
    replicate. NOTE: this is NOT a draw-level resample — it drops the
    replicate's draw from the SBC rank distribution, which biases ranks
    when violations are non-uniform across the posterior. Use only when
    the bounds are soft.

  - `truncate = TRUE`: clamp out-of-bounds response values to the
    nearest bound. This also biases the rank distribution (toward the
    bounds) but keeps the replicate. Neither option implements the
    plan's deterministic draw-resampling; both are documented honestly
    here. A future version may add `on_violation = "resample"`.

- truncate:

  Logical; if `TRUE` and bounds are set, clamp out-of-bounds response
  values instead of producing NAs. Default `FALSE`.

## Value

A generator function `(data_spec, task_ctx) -> data_bundle`.

## Details

Unlike
[`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md)
(which samples from the prior), IFS samples from a preconditioning
posterior, concentrating the truth draw in a region of high posterior
mass. This is the canonical SBC generator for models with diffuse or
improper priors.
