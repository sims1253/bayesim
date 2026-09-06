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
  name = "rmse_loo",
  needs = c("loo", "epred"),
  required = FALSE,
  summary_type = "mean",
  schema = list(value = list(role = "estimate", aggregation = "mean", mcse = "sd"), elpd
    = list(role = "estimate", aggregation = "mean", mcse = "sd"), pareto_k_max =
    list(role = "diagnostic", aggregation = "mean", mcse = "sd"))
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
#>  @ needs       : chr [1:2] "loo" "epred"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ schema      :List of 3
#>  .. $ value       :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ elpd        :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ pareto_k_max:List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
```
