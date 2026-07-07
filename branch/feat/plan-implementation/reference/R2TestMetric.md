# R-squared Test-Set Metric

Coefficient of determination of posterior-mean predictions on the
held-out test set. Returns NA when no test set.

Constructor for R2TestMetric.

## Usage

``` r
R2TestMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
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
```
