# Test-set RMSE compatibility constructor

`rmse_test_metric()` is a naming-compatible wrapper around
[`pred_rmse_metric()`](https://sims1253.github.io/bayesim/reference/RmseMetric.md).
Both compute root-mean-square error of posterior-mean predictions on the
held-out test set; use
[`pred_rmse_metric()`](https://sims1253.github.io/bayesim/reference/RmseMetric.md)
in new code.

## Usage

``` r
rmse_test_metric(name = "rmse_test")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "rmse_test".

## Value

An `RmseMetric` object.

## Examples

``` r
rmse_test_metric()
#> <bayesim::RmseMetric>
#>  @ name        : chr "rmse_test"
#>  @ needs       : chr "predictions"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
