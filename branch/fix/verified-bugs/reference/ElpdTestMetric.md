# ELPD Test-Set Metric

Expected log-predictive density on a held-out test set, estimated by
log-sum-exp over posterior draws of `context$log_lik`. Returns NA when
no test set.

Constructor for ElpdTestMetric.

## Usage

``` r
ElpdTestMetric(
  name = "elpd_test",
  needs = "log_lik",
  required = FALSE,
  summary_type = "mean",
  schema = list(value = list(role = "estimate", aggregation = "mean", mcse = "sd"), n_obs
    = list(role = "count", aggregation = "none", mcse = "none"))
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
#>  @ schema      :List of 2
#>  .. $ value:List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ n_obs:List of 3
#>  ..  ..$ role       : chr "count"
#>  ..  ..$ aggregation: chr "none"
#>  ..  ..$ mcse       : chr "none"
```
