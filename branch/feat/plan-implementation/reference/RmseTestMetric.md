# RMSE Test-Set Metric

Root-mean-square error of posterior-mean predictions on the held-out
test set. Returns NA when no test set.

Constructor for RmseTestMetric.

## Usage

``` r
RmseTestMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)

rmse_test_metric(name = "rmse_test")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "rmse_test".

## Value

A `RmseTestMetric` object.

A `RmseTestMetric` object.

## Examples

``` r
rmse_test_metric()
#> <bayesim::RmseTestMetric>
#>  @ name        : chr "rmse_test"
#>  @ needs       : chr "predictions"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
