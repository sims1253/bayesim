# Execute all tasks in grid

Execute all tasks in grid

## Usage

``` r
execute_tasks(
  task_grid,
  config,
  config_spec,
  fitter,
  metrics,
  retain,
  max_errors,
  progress,
  result_path = NULL,
  config_fingerprint = NULL,
  checkpoint_every = 50L
)
```

## Arguments

- task_grid:

  A task grid tibble from create_task_grid()

- config:

  A SimulationConfig S7 object

- config_spec:

  Plain list config spec for worker transport

- fitter:

  S7 Fitter object

- metrics:

  List of Metric objects

- retain:

  Character vector of what to retain

- max_errors:

  Maximum number of errors before stopping

- progress:

  Logical; if TRUE, show progress bar

- result_path:

  Character; path for checkpoint storage (optional)

- config_fingerprint:

  Character; configuration fingerprint for validation

- checkpoint_every:

  Integer; write checkpoint every N completed tasks. B4: also bounds the
  number of task results held in memory at once.

## Value

A list with task_results and task_grid
