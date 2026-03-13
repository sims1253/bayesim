# bayesim

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![R-CMD-check](https://github.com/sims1253/bayesim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/R-CMD-check.yaml)
[![Tests](https://github.com/sims1253/bayesim/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/test-coverage.yaml)
[![Codecov test coverage](https://codecov.io/gh/sims1253/bayesim/graph/badge.svg)](https://app.codecov.io/gh/sims1253/bayesim)
[![GH-Pages](https://github.com/sims1253/bayesim/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

bayesim provides a simulation framework for reproducible Bayesian modeling studies. It handles task planning, checkpointing, and memory-bounded execution so you can focus on your research questions.

## Installation

You can install the development version of bayesim from GitHub:

```r
# install.packages("pak")
pak::pak("sims1253/bayesim")
```

## Example

The basic workflow has three steps:

1. Create a simulation config with `simulation_config()`
2. Run it with `run_simulation()`
3. Resume interrupted runs with `run_simulation(..., resume = "auto")` or `resume_simulation()`

```r
library(bayesim)

# Define a data generator
data_gen <- function(data_spec, seed, task_ctx) {
  n <- data_spec$n
  x <- rnorm(n)
  y <- data_spec$intercept + data_spec$slope * x + rnorm(n, sd = data_spec$sigma)

  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(
      intercept = data_spec$intercept,
      slope = data_spec$slope,
      sigma = data_spec$sigma
    ),
    vars_of_interest = c("intercept", "slope", "sigma"),
    references = c(intercept = 0, slope = 0, sigma = 1),
    meta = list()
  )
}

# Create the config
config <- simulation_config(
  data_grid = data.frame(
    n = c(100, 500),
    intercept = 1,
    slope = 2,
    sigma = 1
  ),
  fit_grid = data.frame(model = "baseline"),
  data_generator = data_gen,
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 10L,
  seed = 42L
)

# Run the simulation
result <- run_simulation(config, progress = FALSE)
head(result$summary)
```

The engine restores the task RNG state before each call, so repeated, resumed, and parallel runs produce identical results.

## Features

- **Deterministic task planning**: A single study seed determines all task seeds
- **Checkpoint and resume**: Long-running studies can resume after interruption
- **Memory-bounded execution**: `chunk_size` controls how many results stay in memory
- **Extensible design**: S7 classes for fitters and metrics
- **Explicit metrics**: Pass Metric objects instead of string names

## Checkpointing

Set `result_path` and `checkpoint_every` to make runs resumable:

```r
config <- simulation_config(
  data_grid = data.frame(n = c(100, 500)),
  fit_grid = data.frame(model = "baseline"),
  data_generator = data_gen,
  fitter = MockFitter(),
  metrics = list(rmse_metric()),
  n_replicates = 100L,
  seed = 42L,
  result_path = "results/demo-study",
  checkpoint_every = 25L,
  chunk_size = 25L
)

run_simulation(config, resume = "auto")
```

Use `resume_simulation("results/demo-study")` to resume from an existing checkpoint.

## Fitters

bayesim includes:

- `MockFitter()` for testing and examples
- `BrmsFitter()` for brms workflows (default backend: "cmdstanr")

Custom fitters should subclass the `Fitter` S7 class.

## Documentation

See the vignettes for detailed guides:

- `vignette("getting-started")`
- `vignette("simulation-study")`
- `vignette("reproducibility")`
- `vignette("memory-management")`
- `vignette("custom-fitters")`
- `vignette("case-studies")`

## Getting help

If you encounter a bug or have a feature request, please file an issue on [GitHub](https://github.com/sims1253/bayesim/issues).
