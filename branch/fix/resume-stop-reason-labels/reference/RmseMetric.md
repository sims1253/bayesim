# RMSE Metric

Root Mean Square Error between predictions and true values.

Constructor for RmseMetric.

## Usage

``` r
RmseMetric(
  name = "rmse",
  needs = "predictions",
  required = FALSE,
  summary_type = "mean",
  schema = list(value = list(role = "estimate", aggregation = "mean", mcse = "sd"), n_obs
    = list(role = "count", aggregation = "none", mcse = "none"))
)

pred_rmse_metric(name = "rmse")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "rmse".

- needs:

  Character vector of required capabilities. Defaults to "predictions".

- required:

  Logical indicating if metric failure causes task failure. Defaults to
  FALSE.

## Value

An RmseMetric object.

## Examples

``` r
pred_rmse_metric()
#> <bayesim::RmseMetric>
#>  @ name        : chr "rmse"
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
