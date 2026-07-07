# Sampler Diagnostics Metric

Surfaces sampler-level diagnostics (divergences, max treedepth) from
`fit_result$diagnostics`.

Constructor for SamplerDiagnosticsMetric.

## Usage

``` r
SamplerDiagnosticsMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
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
```
