# Run simulation for a single data generation configuration

This function generates multiple datasets according to the data
generation configuration and fits all specified models to each dataset.

## Usage

``` r
dataset_conf_sim(
  data_gen_conf,
  fit_confs,
  prefits,
  seed = NULL,
  result_path = NULL,
  stan_pars,
  ncores,
  cluster_type,
  debug,
  global_seed,
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

- seed:

  Random seed for reproducible data generation

- result_path:

  Path where results and debug files should be saved

- stan_pars:

  A list containing Stan/brms fitting parameters (backend, chains, iter,
  warmup, etc.)

- ncores:

  Number of cores to use for parallel processing

- cluster_type:

  Type of cluster for parallel processing (e.g., "PSOCK")

- debug:

  Logical; if TRUE, saves intermediate results for debugging

- global_seed:

  Global seed for the entire simulation run

- ...:

  Additional arguments passed to dataset_sim and metric calculation
  functions

## Value

A tibble containing simulation results for all datasets and models

## Examples

``` r
if (FALSE) { # \dontrun{
# This function is part of the legacy simulation pipeline.
# New workflows should use simulation_config() + run_simulation().
} # }
```
