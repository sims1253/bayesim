# Run a Simulation Study

Executes a complete simulation study with deterministic reproducibility.

## Usage

``` r
run_simulation(config, resume = c("auto", "never", "must"), progress = TRUE)
```

## Arguments

- config:

  A SimulationConfig S7 object

- resume:

  Character strategy: "auto", "never", or "must"

- progress:

  Logical; if TRUE, show progress bar

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
