# Derive a stable, metric-specific RNG seed

Isolates stochastic metrics from one another: adding or reordering a
metric cannot change another metric's random draws within the same task.

## Usage

``` r
derive_metric_seed(task_seed, metric_name)
```

## Arguments

- task_seed:

  Scalar task seed.

- metric_name:

  Character metric name.

## Value

A positive scalar integer seed.
