# ELPD-LOO Metric

Expected log-predictive-density via PSIS-LOO from `context$loo`.

Constructor for ElpdLooMetric.

## Usage

``` r
ElpdLooMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)

elpd_loo_metric(name = "elpd_loo")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "elpd_loo".

## Value

An `ElpdLooMetric` object.

An `ElpdLooMetric` object.

## Examples

``` r
elpd_loo_metric()
#> <bayesim::ElpdLooMetric>
#>  @ name        : chr "elpd_loo"
#>  @ needs       : chr "loo"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
