# Get Task Count Summary

Returns a summary of task counts by status.

## Usage

``` r
get_task_count_summary(task_grid)
```

## Arguments

- task_grid:

  A task grid tibble.

## Value

Named integer vector with counts for each status.

## Examples

``` r
if (FALSE) { # \dontrun{
summary <- get_task_count_summary(task_grid)
# summary["pending"]  # Number of pending tasks
# summary["success"]  # Number of successful tasks
} # }
```
