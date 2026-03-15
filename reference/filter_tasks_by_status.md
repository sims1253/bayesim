# Filter Tasks by Status

Filters a task grid to include only tasks with specified statuses.

## Usage

``` r
filter_tasks_by_status(task_grid, status)
```

## Arguments

- task_grid:

  A task grid tibble.

- status:

  Character vector of statuses to include. Valid statuses: "pending",
  "success", "failed", "skipped".

## Value

A filtered task grid tibble.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get all completed or skipped tasks
done <- filter_tasks_by_status(task_grid, c("success", "skipped"))
} # }
```
