# Format task ID from indices

Creates a standardized task ID string from data, fit, and replication
indices.

## Usage

``` r
format_task_id(data_idx, fit_idx, rep_idx)
```

## Arguments

- data_idx:

  Integer. Data index (1-999)

- fit_idx:

  Integer. Fit index (1-999)

- rep_idx:

  Integer. Replication index (1-99999)

## Value

Character string in format "dXXX_fXXX_rXXXXX"

## Examples

``` r
format_task_id(1, 2, 100)
#> [1] "d001_f002_r00100"
# Returns: "d001_f002_r00100"
```
