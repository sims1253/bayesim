# Execute a task with safe error handling

Wraps task execution in a tryCatch to ensure recoverable errors are
converted to task results. Fatal errors (bayesim_config_error,
bayesim_contract_error, etc.) are re-thrown to stop the simulation run.

## Usage

``` r
run_task_safe(
  task,
  config_spec,
  fitter,
  metrics,
  retain = c("metrics", "diagnostics")
)
```

## Arguments

- task:

  A task specification list from the task grid.

- config_spec:

  Plain list config spec (for worker transport)

- fitter:

  S7 Fitter object

- metrics:

  List of Metric objects

- retain:

  Character vector of what to retain

## Value

A bayesim_task_result S3 object. If a recoverable error occurs, returns
a failed task result with error information. Fatal errors are re-thrown
and will stop the simulation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run a task safely with error handling
result <- run_task_safe(
  task = task_spec,
  config_spec = config_spec,
  fitter = my_fitter,
  metrics = list(rmse_metric),
  retain = c("metrics", "diagnostics")
)
} # }
```
