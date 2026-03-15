# Run simulation for a single dataset configuration

This function generates a single dataset and fits all specified models
to it, calculating metrics for each model fit.

## Usage

``` r
dataset_sim(
  data_gen_conf,
  fit_confs,
  prefits,
  stan_pars,
  seed,
  debug,
  result_path,
  data_generator_fn = NULL,
  ...
)
```

## Arguments

- data_gen_conf:

  A list containing data generation configuration parameters

- fit_confs:

  A data.frame containing model fitting configurations

- prefits:

  A list of precompiled brms model objects

- stan_pars:

  A list containing Stan/brms fitting parameters

- seed:

  Random seed for reproducible data generation and model fitting

- debug:

  Logical; if TRUE, saves intermediate results for debugging

- result_path:

  Path where debug files should be saved

- data_generator_fn:

  A function with signature `(config, seed) -> list` that generates
  data. The returned list must contain: `dataset`, `sampling_loops`,
  `bad_samples`, `testing_data`, and `true_parameters`.

- ...:

  Additional arguments passed to metric calculation functions

## Value

A tibble containing simulation results for all models

## Examples

``` r
if (FALSE) { # \dontrun{
# This function is part of the legacy simulation pipeline.
# New workflows should use simulation_config() + run_simulation().
} # }
```
