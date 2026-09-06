# Enrich summary with grid columns

Adds data_grid columns (prefixed with "data\_"), fit_grid columns
(prefixed with "fit\_"), and rep_idx to the summary tibble.

## Usage

``` r
enrich_summary_with_grid_columns(summary, task_grid, data_grid, fit_grid)
```

## Arguments

- summary:

  A tibble with task results

- task_grid:

  The task grid tibble with data_idx, fit_idx, rep_idx columns

- data_grid:

  A data frame with data configuration rows

- fit_grid:

  A data frame with fit configuration rows

## Value

The summary tibble with additional columns
