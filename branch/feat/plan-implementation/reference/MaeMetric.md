# MAE Metric

Mean absolute error between predictions and the observed response on the
test set (or training set when no test set).

Constructor for MaeMetric.

## Usage

``` r
MaeMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
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
```
