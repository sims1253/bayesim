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
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
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
#>  @ needs       : chr "loo"
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
