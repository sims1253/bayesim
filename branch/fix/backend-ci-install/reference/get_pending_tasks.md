# Get Pending Tasks

Convenience function to retrieve all tasks that have not yet been
processed.

## Usage

``` r
get_pending_tasks(task_grid)
```

## Arguments

- task_grid:

  A task grid tibble.

## Value

A task grid tibble containing only tasks with status "pending".

## Examples

``` r
if (FALSE) { # \dontrun{
pending <- get_pending_tasks(task_grid)
nrow(pending)  # Number of remaining tasks
} # }
```
