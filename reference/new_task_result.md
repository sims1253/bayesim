# Create a new bayesim_task_result object

Constructs a result object from a single simulation task execution.

## Usage

``` r
new_task_result(
  task_id = character(1),
  status = "success",
  metrics = NULL,
  diagnostics = NULL,
  timing = list(total = 0),
  error = NULL,
  warnings = character()
)
```

## Arguments

- task_id:

  Character scalar identifying the task

- status:

  Character scalar: one of "success", "failed", or "skipped"

- metrics:

  Named list of computed metrics (NULL if task failed or skipped)

- diagnostics:

  Named list of diagnostic values (NULL if task failed)

- timing:

  List containing timing information, must include `total`

- error:

  NULL, or a list with `error_class` and `error_message` if failed

- warnings:

  Character vector of warning messages

## Value

A validated `bayesim_task_result` object

## Details

The bayesim_task_result class captures the outcome of a single
simulation task, including its computed metrics, any diagnostics, timing
information, and errors or warnings.

Validation rules:

- `status` must be one of "success", "failed", or "skipped"

- If `status` is "success", `metrics` must not be NULL

- If `status` is "failed", `error` must not be NULL

- `timing$total` must be non-negative

## Examples

``` r
# Successful task
result <- new_task_result(
  task_id = "task_001",
  status = "success",
  metrics = list(rmse = 0.05, bias = 0.01),
  diagnostics = list(n_eff = 500, rhat = 1.01),
  timing = list(total = 5.2)
)

# Failed task
result <- new_task_result(
  task_id = "task_002",
  status = "failed",
  error = list(error_class = "convergence_error", error_message = "R-hat > 1.1"),
  timing = list(total = 2.0)
)
```
