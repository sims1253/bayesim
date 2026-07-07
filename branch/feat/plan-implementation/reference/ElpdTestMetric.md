# ELPD Test-Set Metric

Expected log-predictive density on a held-out test set, estimated by
log-sum-exp over posterior draws of `context$log_lik`. Returns NA when
no test set.

Constructor for ElpdTestMetric.

## Usage

``` r
ElpdTestMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)

elpd_test_metric(name = "elpd_test")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "elpd_test".

## Value

An `ElpdTestMetric` object.

An `ElpdTestMetric` object.

## Examples

``` r
elpd_test_metric()
#> <bayesim::ElpdTestMetric>
#>  @ name        : chr "elpd_test"
#>  @ needs       : chr "log_lik"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
