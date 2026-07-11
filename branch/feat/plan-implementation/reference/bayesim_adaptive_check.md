# Should the run stop early based on the adaptive-stopping policy? (I3)

Builds a quick summary from `task_results_so_far` and calls
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)
for `stop_on$estimand`. Returns TRUE when the MCSE of `stop_on$measure`
is finite and below `stop_on$target_mcse` in every condition cell.
Wrapped in tryCatch by the caller: any failure (e.g. no truth columns
yet) =\> FALSE.

## Usage

``` r
bayesim_adaptive_check(task_results_so_far, task_grid, stop_on, config)
```

## Arguments

- task_results_so_far:

  List of in-memory `bayesim_task_result`.

- task_grid:

  Current task grid (statuses updated).

- stop_on:

  A validated `stop_on` list (see
  [`validate_stop_on()`](https://sims1253.github.io/bayesim/reference/validate_stop_on.md)).

- config:

  A SimulationConfig.

## Value

TRUE/FALSE. Never throws.
