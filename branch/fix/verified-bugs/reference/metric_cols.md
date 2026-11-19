# Select flattened metric columns

Selects columns belonging to one metric in bayesim's flattened naming
scheme, `<metric>__<field>` or `<metric>__<field>__<param>`.

## Usage

``` r
metric_cols(x, metric, fields = NULL, as = c("names", "long"))
```

## Arguments

- x:

  A bayesim simulation result (with a `summary` component), or a
  per-task summary data frame.

- metric:

  Required character scalar naming the metric prefix.

- fields:

  Optional character vector restricting the field segment.

- as:

  Output form: `"names"` returns a named character vector suitable for
  [`dplyr::all_of()`](https://tidyselect.r-lib.org/reference/all_of.html);
  `"long"` returns one row per task and metric value.

## Value

For `as = "names"`, a named character vector whose values are full
column names and whose names are suffixes following the metric prefix.
For `as = "long"`, a tibble with optional `task_id`, `field`, `param`,
and `value` columns.

## Examples

``` r
fitter <- LinearRegressionFitter()
summary <- data.frame(
  task_id = "task_1",
  posterior_summary__mean__x = 1.2,
  posterior_summary__sd__x = 0.3,
  check.names = FALSE
)
metric_cols(summary, "posterior_summary")
#>                      mean__x                        sd__x 
#> "posterior_summary__mean__x"   "posterior_summary__sd__x" 
metric_cols(summary, "posterior_summary", fields = "mean", as = "long")
#> # A tibble: 1 × 4
#>   task_id field param value
#>   <chr>   <chr> <chr> <dbl>
#> 1 task_1  mean  x       1.2
```
