# MSE Metric

Mean squared error between predictions and the observed response on the
test set (or training set when no test set).

Constructor for MseMetric.

## Usage

``` r
MseMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)

pred_mse_metric(name = "mse")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "mse".

## Value

An `MseMetric` object.

An `MseMetric` object.

## Examples

``` r
pred_mse_metric()
#> <bayesim::MseMetric>
#>  @ name        : chr "mse"
#>  @ needs       : chr "predictions"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
