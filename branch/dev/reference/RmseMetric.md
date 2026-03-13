# RMSE Metric

Root Mean Square Error between predictions and true values.

`rmse_metric()` is a constructor function that creates an RmseMetric
instance with appropriate defaults. Use this constructor to work around
S7's property default inheritance issue.

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
rmse_metric(name = "my_rmse")
#> <bayesim::RmseMetric>
#>  @ name    : chr "my_rmse"
#>  @ needs   : chr "predictions"
#>  @ required: logi FALSE
```
