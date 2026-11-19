# Convert task results to data frame

Converts a list of `bayesim_task_result` objects to a data frame
suitable for checkpoint storage.

## Usage

``` r
results_to_dataframe(task_results)
```

## Arguments

- task_results:

  A list of `bayesim_task_result` objects.

## Value

A data frame with one row per task, containing:

- `task_id` - Task identifier

- `status` - Task status ("success", "failed", "skipped")

- Columns from `metrics` (if present)

- Columns from `diagnostics` (if present)

## Details

This function flattens the nested structure of task results into a
tabular format. Metrics and diagnostics are extracted and added as
columns. NULL task results are skipped.
