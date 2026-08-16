# Build a quick summary tibble from in-memory task results (I3)

Internal helper for adaptive stopping: flattens the non-NULL entries of
`task_results` to a summary data.frame via
[`results_to_dataframe()`](https://sims1253.github.io/bayesim/reference/results_to_dataframe.md)
and enriches it with data_grid/fit_grid/rep_idx columns (matching what
[`build_simulation_result()`](https://sims1253.github.io/bayesim/reference/build_simulation_result.md)
produces). Used by the internal adaptive evaluator so it can call
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)
mid-run.

## Usage

``` r
bayesim_adaptive_summary(task_results, task_grid, config)
```

## Arguments

- task_results:

  List of `bayesim_task_result` (possibly with NULLs).

- task_grid:

  The task grid tibble (with up-to-date statuses).

- config:

  A SimulationConfig.

## Value

A data.frame summary, or NULL if no non-NULL results.
