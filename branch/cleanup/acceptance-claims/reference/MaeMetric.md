# MAE Metric

Mean absolute error between predictions and the observed response on the
test set. Returns NA when no test set is present (no training-set
fallback).

Constructor for MaeMetric.

## Usage

``` r
MaeMetric(
  name = "mae",
  needs = "predictions",
  required = FALSE,
  summary_type = "mean",
  schema = list(value = list(role = "estimate", aggregation = "mean", mcse = "sd"), n_obs
    = list(role = "count", aggregation = "none", mcse = "none"))
)

pred_mae_metric(name = "mae")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "mae".

## Value

A `MaeMetric` object.

A `MaeMetric` object.

## Examples

``` r
pred_mae_metric()
#> <bayesim::MaeMetric>
#>  @ name        : chr "mae"
#>  @ needs       : chr "predictions"
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
