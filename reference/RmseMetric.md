# RMSE Metric

Root Mean Square Error between predictions and true values.

Constructor for RmseMetric.

## Usage

``` r
RmseMetric(name = character(0), needs = character(0), required = FALSE)

rmse_metric(name = "rmse")
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
rmse_metric()
#> <bayesim::RmseMetric>
#>  @ name    : chr "rmse"
#>  @ needs   : chr "predictions"
#>  @ required: logi FALSE
```
