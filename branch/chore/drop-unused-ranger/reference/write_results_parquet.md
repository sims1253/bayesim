# Write a simulation summary to parquet

Writes `results_df` (a data frame/tibble summary produced by
[`build_simulation_result()`](https://sims1253.github.io/bayesim/reference/build_simulation_result.md))
to `path` as a parquet file using `nanoparquet`. The data frame is
coerced to a plain data.frame first so that tibble/vctrs columns parquet
cannot represent are flattened.

## Usage

``` r
write_results_parquet(results_df, path)
```

## Arguments

- results_df:

  A data.frame or tibble summary to write.

- path:

  Character scalar. Destination parquet file path.

## Value

Invisible `path` (invisibly), the path written to.

## Details

Requires the suggested `nanoparquet` package; if it is not installed,
[`rlang::check_installed()`](https://rlang.r-lib.org/reference/is_installed.html)
prompts (or errors) accordingly.
