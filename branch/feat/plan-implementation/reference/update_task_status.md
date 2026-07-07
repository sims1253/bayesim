# Update Task Status in Grid

Updates the status of a single task in the task grid. This function is
designed for functional use (returns modified copy).

## Usage

``` r
update_task_status(task_grid, task_id, status)
```

## Arguments

- task_grid:

  A task grid tibble.

- task_id:

  The task ID to update.

- status:

  New status value. Valid values: "success", "failed", "skipped".

## Value

A modified task grid tibble with the updated status.

## Examples

``` r
if (FALSE) { # \dontrun{
task_grid <- update_task_status(task_grid, "d001_f001_r00001", "success")
} # }
```
