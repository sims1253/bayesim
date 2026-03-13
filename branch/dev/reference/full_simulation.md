# Run a full simulation study

Legacy interface for older bayesim simulation studies.

## Usage

``` r
full_simulation(
  data_gen_confs,
  fit_confs,
  ncores_simulation = 1,
  cluster_type = "PSOCK",
  stan_pars,
  seed = NULL,
  result_path = NULL,
  debug = FALSE,
  ...
)
```

## Arguments

- data_gen_confs:

  A data.frame where each row specifies a data generation configuration

- fit_confs:

  A data.frame where each row specifies a model fitting configuration

- ncores_simulation:

  Number of cores to use for parallel dataset simulation

- cluster_type:

  Type of cluster for parallel processing (e.g., "PSOCK", "FORK")

- stan_pars:

  A list containing Stan/brms fitting parameters:

  - backend: Stan backend ("cmdstanr" or "rstan")

  - chains: Number of MCMC chains

  - iter: Number of iterations per chain

  - warmup: Number of warmup iterations

  - init: Initial values for MCMC

  - cmdstan_path: Path to CmdStan (if using cmdstanr backend)

  - cmdstan_write_path: Path for CmdStan output files

- seed:

  Random seed for reproducible simulation

- result_path:

  Path where results should be saved (as RDS files)

- debug:

  Logical; if TRUE, enables verbose output and saves intermediate
  results

- ...:

  Additional arguments passed to dataset_conf_sim and metric calculation
  functions

## Value

A tibble containing all simulation results

## Details

New code should prefer
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md),
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md),
and
[`resume_simulation()`](https://sims1253.github.io/bayesim/reference/resume_simulation.md).
This wrapper is retained for compatibility with the pre-rewrite API.

## Examples

``` r
if (FALSE) { # \dontrun{
# Define data generation configurations
data_gen_confs <- data.frame(
  id = c("config1", "config2"),
  data_N = c(100, 100),
  dataset_N = c(10, 10)
)

# Define model fitting configurations
fit_confs <- data.frame(
  fit_family = c("gaussian", "student_t"),
  fit_link = c("identity", "identity"),
  formula = c("y ~ x", "y ~ x")
)

# Set Stan parameters
stan_pars <- list(
  backend = "cmdstanr",
  chains = 4,
  iter = 2000,
  warmup = 1000,
  init = 0.1
)

# Run simulation
results <- full_simulation(
  data_gen_confs = data_gen_confs,
  fit_confs = fit_confs,
  stan_pars = stan_pars,
  ncores_simulation = 4,
  seed = 12345
)
} # }
```
