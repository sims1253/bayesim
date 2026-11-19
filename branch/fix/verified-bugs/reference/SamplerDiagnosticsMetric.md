# Sampler Diagnostics Metric

Surfaces the full set of sampler and convergence diagnostics from
`fit_result$diagnostics`: max R-hat, min bulk ESS, min tail ESS,
divergence count, and max-treedepth hits. Returns NA fields when the
diagnostics list is absent.

Constructor for SamplerDiagnosticsMetric.

## Usage

``` r
SamplerDiagnosticsMetric(
  name = "sampler_diagnostics",
  needs = character(0),
  required = FALSE,
  summary_type = "mean",
  schema = list(rhat_max = list(role = "diagnostic", aggregation = "mean", mcse = "sd"),
    ess_bulk_min = list(role = "diagnostic", aggregation = "mean", mcse = "sd"),
    ess_tail_min = list(role = "diagnostic", aggregation = "mean", mcse = "sd"),
    divergent = list(role = "count", aggregation = "none", mcse = "none"), max_treedepth
    = list(role = "diagnostic", aggregation = "mean", mcse = "sd"))
)

sampler_diagnostics_metric(name = "sampler_diagnostics")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "sampler_diagnostics".

## Value

A `SamplerDiagnosticsMetric` object.

A `SamplerDiagnosticsMetric` object.

## Examples

``` r
sampler_diagnostics_metric()
#> <bayesim::SamplerDiagnosticsMetric>
#>  @ name        : chr "sampler_diagnostics"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ schema      :List of 5
#>  .. $ rhat_max     :List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ ess_bulk_min :List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ ess_tail_min :List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ divergent    :List of 3
#>  ..  ..$ role       : chr "count"
#>  ..  ..$ aggregation: chr "none"
#>  ..  ..$ mcse       : chr "none"
#>  .. $ max_treedepth:List of 3
#>  ..  ..$ role       : chr "diagnostic"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
```
