# Create a new bayesim_simulation_result object

Constructs a result object from a complete simulation run.

## Usage

``` r
new_simulation_result(
  config_fingerprint = character(1),
  task_results = list(),
  task_grid = NULL,
  summary = NULL,
  timing = list(total = 0),
  errors = NULL,
  checkpoint_path = NULL,
  metric_summary_types = NULL,
  metric_field_metadata = NULL
)
```

## Arguments

- config_fingerprint:

  Character hash uniquely identifying the configuration

- task_results:

  List of `bayesim_task_result` objects

- task_grid:

  Tibble with task grid information (task_id, data_idx, fit_idx,
  rep_idx, status)

- summary:

  Tibble with one row per task, columns for metrics and diagnostics

- timing:

  List containing total, by_phase, and other timing breakdowns

- errors:

  Tibble of failed tasks with error details

- checkpoint_path:

  Path to the final checkpoint file, or NULL

## Value

A validated `bayesim_simulation_result` object

## Details

The bayesim_simulation_result class encapsulates the complete results of
a simulation study, including all individual task results, an aggregated
summary, timing breakdowns, error information, and checkpoint location.

Validation rules:

- `config_fingerprint` must be a scalar character

- `task_results` must be a list where all elements are
  `bayesim_task_result`

- `summary` and `errors` must be data frames

- `timing$total` must be non-negative

- `checkpoint_path` must be NULL or a scalar character

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a simulation result
task1 <- new_task_result(
  task_id = "task_001",
  status = "success",
  metrics = list(rmse = 0.05),
  timing = list(total = 5.0)
)

result <- new_simulation_result(
  config_fingerprint = "abc123",
  task_results = list(task1),
  task_grid = tibble::tibble(
    task_id = "task_001",
    data_idx = 1L,
    fit_idx = 1L,
    rep_idx = 1L,
    status = "success"
  ),
  summary = tibble::tibble(task_id = "task_001", rmse = 0.05),
  timing = list(total = 10.0, by_phase = list(setup = 2.0, fit = 8.0)),
  errors = tibble::tibble(task_id = character(), error_message = character()),
  checkpoint_path = "/path/to/checkpoint.rds"
)
} # }
```
