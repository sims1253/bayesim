# Convergence Metric

Summarizes convergence diagnostics (max R-hat, min ESS, divergences)
from `fit_result$diagnostics`. Returns NA when absent.

Constructor for ConvergenceMetric.

## Usage

``` r
ConvergenceMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)

convergence_metric(name = "convergence")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "convergence".

## Value

A `ConvergenceMetric` object.

A `ConvergenceMetric` object.

## Examples

``` r
convergence_metric()
#> <bayesim::ConvergenceMetric>
#>  @ name        : chr "convergence"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
```
