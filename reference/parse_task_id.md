# Parse Task ID Components

Extracts the data_idx, fit_idx, and rep_idx from a task ID string.

## Usage

``` r
parse_task_id(task_id)
```

## Arguments

- task_id:

  Character task ID in format "dXXX_fXXX_rXXXXX".

## Value

Named integer vector with elements data_idx, fit_idx, rep_idx, or NULL
if the task_id format is invalid.

## Examples

``` r
if (FALSE) { # \dontrun{
components <- parse_task_id("d002_f003_r00042")
# components$data_idx == 2
# components$fit_idx == 3
# components$rep_idx == 42
} # }
```
