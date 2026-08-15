# Getting Started with bayesim

## Introduction

bayesim is a simulation framework for Bayesian modeling studies. It
provides:

- **Reproducible execution** across sequential, parallel, and resumed
  runs
- **Checkpoint/resume** capabilities for long-running simulations
- **Memory-bounded execution** with configurable artifact retention
- **Extensible interfaces** for custom fitters and metrics

## Basic Usage

> **Note:** This vignette uses
> [`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
> which requires the brms package. For quick testing without brms,
> replace
> [`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
> with
> [`MockFitter()`](https://sims1253.github.io/bayesim/reference/MockFitter.md).

### Setting Up a Simulation

``` r

library(bayesim)
```

#### Define a Data Generator

The data generator creates synthetic datasets for your simulation. It
must have the signature `(data_spec, seed, task_ctx)`:

``` r

my_data_generator <- function(data_spec, seed, task_ctx) {
  # Note: seed is a scalar task seed.
  # The simulation engine also restores the task RNG stream before each call,
  # so repeated full, resumed, and parallel runs stay aligned.
  n <- data_spec$n
  
  # Generate synthetic data
  x <- rnorm(n)
  y <- data_spec$intercept + data_spec$slope * x + rnorm(n, sd = data_spec$sigma)
  
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(intercept = data_spec$intercept, slope = data_spec$slope, sigma = data_spec$sigma),
    vars_of_interest = c("intercept", "slope", "sigma"),
    references = c(intercept = 0, slope = 0, sigma = 1),
    meta = list()
  )
}
```

#### Create a Configuration

Create a simulation configuration with your data grid, fit grid, and
other parameters. Note that metrics should be passed as a list of Metric
objects. Built-in metric constructors like
[`rmse_metric()`](https://sims1253.github.io/bayesim/reference/RmseMetric.md)
have default names, but custom metrics require an explicit `name`
argument:

``` r

# Note: brms must be installed to use BrmsFitter()
config <- simulation_config(
  data_grid = data.frame(
    n = c(100, 500),
    intercept = 1,
    slope = 2,
    sigma = 1
  ),
  fit_grid = data.frame(
    model = "linear"
  ),
  data_generator = my_data_generator,
  fitter = BrmsFitter(),  # Use BrmsFitter for Bayesian model fitting
  metrics = list(
    rmse_metric(),
    bias_metric()
  ),
  n_replicates = 10L,
  seed = 42L
)
```

``` r

# Alternative: Use MockFitter for testing without brms
# Note: MockFitter is for testing the simulation framework only.
# For real Bayesian inference, use BrmsFitter() or a custom fitter.
# See vignette("custom-fitters") for examples.
config <- simulation_config(
  data_grid = data.frame(
    n = c(100, 500),
    intercept = 1,
    slope = 2,
    sigma = 1
  ),
  fit_grid = data.frame(
    model = "linear"
  ),
  data_generator = my_data_generator,
  fitter = MockFitter(),  # Use mock for quick testing
  metrics = list(
    rmse_metric(),
    bias_metric()
  ),
  n_replicates = 10L,
  seed = 42L
)
```

#### Run the Simulation

``` r

result <- run_simulation(config, progress = FALSE)
```

#### Examine Results

``` r

print(result)

# Summary tibble
head(result$summary)
```

The summary tibble includes:

- `task_id`, `status`, `timing_total`: Basic task information
- `data_<colname>`: Columns from your data_grid (e.g., `data_n`,
  `data_sigma`)
- `fit_<colname>`: Columns from your fit_grid (e.g., `fit_model`)
- `<metric>__<field>`: Metric outputs (e.g., `rmse__value`,
  `bias__value`)

## Checkpointing and Resume

For long-running simulations, use checkpointing:

``` r

config <- simulation_config(
  data_grid = my_data_grid,
  fit_grid = my_fit_grid,
  data_generator = my_data_generator,
  fitter = my_fitter,
  metrics = list(rmse_metric()),
  n_replicates = 1000L,
  seed = 42L,
  result_path = "my_simulation",
  checkpoint_every = 50L,
  chunk_size = 50L,
  checkpoint_format = "rds"
)

# Run (can be interrupted and resumed)
result <- run_simulation(config, resume = "auto")
```

If interrupted, resume with:

``` r

result <- resume_simulation("my_simulation")
```

`checkpoint_format = "parquet"` is reserved for a future backend but is
not implemented yet.

## Custom Metrics

Create custom metrics by extending the Metric class:

``` r

MyMetric <- S7::new_class(
  "MyMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "my_metric"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

S7::method(compute, MyMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  # Your metric computation
  list(value = 0.5)
}

# Use directly in simulation_config
config <- simulation_config(
  data_grid = data.frame(n = 100),
  fit_grid = data.frame(model = "test"),
  data_generator = my_data_generator,
  fitter = MockFitter(),
  metrics = list(MyMetric(name = "my_metric")),  # Use your custom metric with explicit name
  n_replicates = 5L,
  seed = 42L
)
```

## Next Steps

- See
  [`vignette("simulation-study")`](https://sims1253.github.io/bayesim/articles/simulation-study.md)
  for a complete example simulation study with analysis and
  visualizations
- See
  [`vignette("custom-fitters")`](https://sims1253.github.io/bayesim/articles/custom-fitters.md)
  for creating custom model fitters
- See
  [`vignette("reproducibility")`](https://sims1253.github.io/bayesim/articles/reproducibility.md)
  for understanding determinism guarantees
- See
  [`vignette("memory-management")`](https://sims1253.github.io/bayesim/articles/memory-management.md)
  for handling large simulations
