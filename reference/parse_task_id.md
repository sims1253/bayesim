# Parse task ID to indices

Parses a task ID string back into its component indices.

Extracts the data_idx, fit_idx, and rep_idx from a task ID string.

## Usage

``` r
parse_task_id(task_id)

parse_task_id(task_id)
```

## Arguments

- task_id:

  Character task ID in format "dXXX_fXXX_rXXXXX".

## Value

Named list with elements `data_idx`, `fit_idx`, and `rep_idx`.

Named integer vector with elements data_idx, fit_idx, rep_idx, or NULL
if the task_id format is invalid.

## Examples

``` r
parse_task_id("d001_f002_r00100")
#> $data_idx
#> [1] 1
#> 
#> $fit_idx
#> [1] 2
#> 
#> $rep_idx
#> [1] 100
#> 
# Returns: list(data_idx = 1, fit_idx = 2, rep_idx = 100)
if (FALSE) { # \dontrun{
components <- parse_task_id("d002_f003_r00042")
# components$data_idx == 2
# components$fit_idx == 3
# components$rep_idx == 42
} # }
```
