# RMSE-LOO Metric

RMSE estimated via PSIS-LOO: the LOO-weighted posterior-mean prediction
is computed with
[`loo::E_loo()`](https://mc-stan.org/loo/reference/E_loo.html) and
compared to the observed response (`sqrt(mean((y - yloo)^2))`). This
matches the construction underlying brms' `loo_predict(type = "mean")`.
The max Pareto-k is returned as a diagnostic. Falls back to NA when the
PSIS object or observed response is unavailable.

Constructor for RmseLooMetric.

## Usage

``` r
RmseLooMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)

rmse_loo_metric(name = "rmse_loo")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "rmse_loo".

## Value

A `RmseLooMetric` object.

A `RmseLooMetric` object.

## Examples

``` r
rmse_loo_metric()
#> <bayesim::RmseLooMetric>
#>  @ name        : chr "rmse_loo"
#>  @ needs       : chr "loo"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
