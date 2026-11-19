# Posterior Summary Metric

Posterior mean, median, standard deviation, and quantile-based credible
intervals for each parameter in `vars_of_interest`.

Constructor for PosteriorSummaryMetric.

## Usage

``` r
PosteriorSummaryMetric(
  name = "posterior_summary",
  needs = character(0),
  required = FALSE,
  summary_type = "mean",
  schema = list(mean = list(role = "estimate", aggregation = "mean", mcse = "sd"), median
    = list(role = "estimate", aggregation = "mean", mcse = "sd"), sd = list(role =
    "estimate", aggregation = "mean", mcse = "sd"), q_lower = list(role = "estimate",
    aggregation = "mean", mcse = "sd"), q_upper = list(role = "estimate", aggregation =
    "mean", mcse = "sd")),
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
#>  @ schema      :List of 5
#>  .. $ mean   :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ median :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ sd     :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ q_lower:List of 4
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  ..  ..$ nominal    : num 0.95
#>  .. $ q_upper:List of 4
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  ..  ..$ nominal    : num 0.95
#>  @ prob        : num 0.95
```
