# ELPD-LOO Metric

Expected log-predictive-density via PSIS-LOO from `context$loo`.

Constructor for ElpdLooMetric.

## Usage

``` r
ElpdLooMetric(
  name = "elpd_loo",
  needs = "loo",
  required = FALSE,
  summary_type = "mean",
  schema = list(elpd = list(role = "estimate", aggregation = "mean", mcse = "sd"), p_loo
    = list(role = "estimate", aggregation = "mean", mcse = "sd"), se = list(role =
    "estimate", aggregation = "mean", mcse = "sd"), pareto_k_max = list(role =
    "diagnostic", aggregation = "mean", mcse = "sd"))
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
#>  @ schema      :List of 4
#>  .. $ elpd        :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ p_loo       :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ se          :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ pareto_k_max:List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
```
