# Read a simulation summary file

Reads a simulation summary written by
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
(or
[`write_results_parquet()`](https://sims1253.github.io/bayesim/reference/write_results_parquet.md)).
Dispatches on file extension: `.parquet` files are read with
[`nanoparquet::read_parquet`](https://nanoparquet.r-lib.org/reference/read_parquet.html)
(requires the suggested `nanoparquet` package); `.rds` files are read
with [`base::readRDS()`](https://rdrr.io/r/base/readRDS.html).

## Usage

``` r
read_summary(path)
```

## Arguments

- path:

  Character scalar. Path to a `.parquet` or `.rds` summary file.

## Value

A data frame of the summary.

## Examples

``` r
if (FALSE) { # \dontrun{
summary_df <- read_summary("results/summary.parquet")
summary_df <- read_summary("results/summary.rds")
} # }
```
