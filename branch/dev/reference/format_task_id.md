# Format task ID from indices (deprecated)

**\[deprecated\]**

Use `make_task_id()` instead, which auto-calculates field widths from
the grid dimensions.

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
#> Warning: `format_task_id()` was deprecated in bayesim 1.1.
#> ℹ Please use `make_task_id()` instead.
#> [1] "d001_f002_r00100"
# Returns: "d001_f002_r00100"
```
