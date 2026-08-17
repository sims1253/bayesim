# R-squared Test-Set Metric

Coefficient of determination of posterior-mean predictions on the
held-out test set. Returns NA when no test set.

Constructor for R2TestMetric.

## Usage

``` r
R2TestMetric(
  name = "r2_test",
  needs = "predictions",
  required = FALSE,
  summary_type = "mean",
  schema = list(value = list(role = "estimate", aggregation = "mean", mcse = "sd"), n_obs
    = list(role = "count", aggregation = "none", mcse = "none"), undefined = list(role =
    "diagnostic", aggregation = "none", mcse = "none"))
)

r2_test_metric(name = "r2_test")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "r2_test".

## Value

A `R2TestMetric` object.

A `R2TestMetric` object.

## Examples

``` r
r2_test_metric()
#> <bayesim::R2TestMetric>
#>  @ name        : chr "r2_test"
#>  @ needs       : chr "predictions"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ schema      :List of 3
#>  .. $ value    :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ n_obs    :List of 3
#>  ..  ..$ role       : chr "count"
#>  ..  ..$ aggregation: chr "none"
#>  ..  ..$ mcse       : chr "none"
#>  .. $ undefined:List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "none"
#>  ..  ..$ mcse       : chr "none"
```
