# Run a Simulation Study

Executes a complete simulation study with deterministic reproducibility.

## Usage

``` r
run_simulation(
  config,
  resume = c("auto", "never", "must"),
  progress = TRUE,
  workers = NULL
)
```

## Arguments

- config:

  A SimulationConfig S7 object

- resume:

  Character strategy: "auto", "never", or "must"

- progress:

  Logical; if TRUE, show progress bar

- workers:

  Positive integer, NULL, or "multisession". When non-NULL,
  `mirai::daemons(workers)` is set up for the run and torn down on exit
  — the simple path for local parallelism. Must be NULL when daemons are
  already set (use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  directly for the advanced/HPC path: remote daemons, TLS, etc.).
  Daemons are set before the model bank ships.

## Value

A bayesim_simulation_result S3 object

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(
  data_grid = data.frame(n = c(100, 500)),
  fit_grid = data.frame(model = "baseline"),
  data_generator = my_data_gen,
  fitter = my_fitter,
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 100L,
  seed = 42L
)
result <- run_simulation(config)
} # }
```
