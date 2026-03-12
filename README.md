# Bayesim

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![R-CMD-check](https://github.com/sims1253/bayesim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/R-CMD-check.yaml)
[![Tests](https://github.com/sims1253/bayesim/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/test-coverage.yaml)
[![Codecov test coverage](https://codecov.io/gh/sims1253/bayesim/graph/badge.svg)](https://app.codecov.io/gh/sims1253/bayesim)
[![GH-Pages](https://github.com/sims1253/bayesim/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/pkgdown.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

`bayesim` is a simulation framework for reproducible Bayesian modeling studies.

The current user-facing workflow is:

1. Build a `SimulationConfig` with `simulation_config()`
2. Run it with `run_simulation()`
3. Resume interrupted runs with `run_simulation(..., resume = "auto" | "must")` or `resume_simulation()`

## Core ideas

- Deterministic task planning from a single study seed
- Checkpoint/resume for long-running studies
- Memory-bounded execution via `chunk_size`
- S7 fitters and metrics for extensibility
- Explicit Metric objects instead of string metric names

## Minimal example

```r
library(bayesim)

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

result <- run_simulation(config, progress = FALSE)
head(result$summary)
```

`data_generator()` receives a scalar `seed` and a `task_ctx`. The engine also restores the task RNG state before each call, so repeated full, resumed, and parallel runs stay aligned.

## Checkpointing and resume

Use `result_path` plus `checkpoint_every` to make runs resumable:

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
resume_simulation("results/demo-study")
```

- `checkpoint_format = "rds"` is the supported checkpoint backend
- `checkpoint_format = "parquet"` is reserved but not implemented yet
- `chunk_size` controls how many task results are retained before a flush/checkpoint cycle
- `max_in_memory` still works as a compatibility alias for `chunk_size`

## Metrics

Pass Metric objects, not strings:

```r
metrics = list(
  rmse_metric(),
  bias_metric(),
  coverage_metric()
)
```

Custom metrics should subclass `Metric` and supply a stable `name`.

## Fitters

`bayesim` ships with `MockFitter()` for fast examples and `BrmsFitter()` for `brms` workflows. The default `brms` backend is `"cmdstanr"`.

## Large studies

For larger simulations:

- set `result_path`
- choose `checkpoint_every`
- set `chunk_size` to cap in-memory task results
- use a minimal retention profile when heavy artifacts are unnecessary

High-cardinality metric payloads are written to artifact files instead of expanding the main summary table indefinitely.

## Vignettes

- `vignette("getting-started")`
- `vignette("simulation-study")`
- `vignette("reproducibility")`
- `vignette("memory-management")`
- `vignette("custom-fitters")`
- `vignette("case-studies")`
