# Fit a single model in a simulation study

This function fits a Bayesian model using a precompiled brms model
object and calculates specified metrics on the fitted model.

## Usage

``` r
fit_sim(
  prefit,
  dataset,
  data_gen_output,
  fit_conf,
  seed,
  debug,
  result_path,
  stan_pars,
  ...
)
```

## Arguments

- prefit:

  A precompiled brms model object from
  [`get_prefit`](https://sims1253.github.io/bayesim/reference/get_prefit.md)

- dataset:

  A data.frame containing the data to fit the model to

- data_gen_output:

  Output from the data generation function containing true parameters

- fit_conf:

  A list or data.frame row containing model configuration (formula,
  family, etc.)

- seed:

  Random seed for reproducible model fitting

- debug:

  Logical; if TRUE, saves intermediate results for debugging

- result_path:

  Path where debug files should be saved

- stan_pars:

  List of Stan/brms parameters (chains, iter, warmup, etc.)

- ...:

  Additional arguments passed to metric calculation functions

## Value

A tibble containing calculated metrics for the fitted model

## Examples

``` r
if (FALSE) { # \dontrun{
# This function is part of the legacy simulation pipeline.
# New workflows should use simulation_config() + run_simulation().
} # }
```
