# R-squared LOO Metric

Variance-explained estimated via PSIS-LOO following brms' `loo_R2()`:
`1 - var_loo(y - yloo) / var_loo(y)`, where `yloo` is the LOO-weighted
posterior expectation
([`loo::E_loo()`](https://mc-stan.org/loo/reference/E_loo.html) mean of
`posterior_epred`) and the variances use the same weighted-expecation
construction as brms (Gelman, Goodrich, Gabry & Vehtari 2018). Falls
back to NA when epred or the PSIS object is unavailable.

Constructor for R2LooMetric.

## Usage

``` r
R2LooMetric(
  name = "r2_loo",
  needs = c("loo", "epred"),
  required = FALSE,
  summary_type = "mean",
  schema = list(value = list(role = "estimate", aggregation = "mean", mcse = "sd"), elpd
    = list(role = "estimate", aggregation = "mean", mcse = "sd"))
)

r2_loo_metric(name = "r2_loo")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "r2_loo".

## Value

An `R2LooMetric` object.

An `R2LooMetric` object.

## Examples

``` r
r2_loo_metric()
#> <bayesim::R2LooMetric>
#>  @ name        : chr "r2_loo"
#>  @ needs       : chr [1:2] "loo" "epred"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ schema      :List of 3
#>  .. $ value    :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ elpd     :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ undefined:List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "none"
#>  ..  ..$ mcse       : chr "none"
```
