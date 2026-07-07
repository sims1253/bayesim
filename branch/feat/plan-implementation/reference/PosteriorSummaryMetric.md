# Posterior Summary Metric

Posterior mean, median, standard deviation, and quantile-based credible
intervals for each parameter in `vars_of_interest`.

Constructor for PosteriorSummaryMetric.

## Usage

``` r
PosteriorSummaryMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean",
  prob = 0.95
)

posterior_summary_metric(name = "posterior_summary", prob = 0.95)
```

## Arguments

- name:

  Character string naming the metric. Defaults to "posterior_summary".

- prob:

  Credible-interval mass. Defaults to 0.95.

## Value

A `PosteriorSummaryMetric` object.

A `PosteriorSummaryMetric` object.

## Examples

``` r
posterior_summary_metric()
#> <bayesim::PosteriorSummaryMetric>
#>  @ name        : chr "posterior_summary"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ prob        : num 0.95
```
