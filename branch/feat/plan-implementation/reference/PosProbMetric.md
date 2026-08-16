# Positivity Probability Metric

For each parameter in `vars_of_interest`, the posterior probability that
the parameter is positive (fraction of draws \> 0).

Constructor for PosProbMetric.

## Usage

``` r
PosProbMetric(
  name = "pos_prob",
  needs = character(0),
  required = FALSE,
  summary_type = "mean",
  schema = list(mean = list(role = "estimate", aggregation = "mean", mcse = "sd"),
    by_param = list(role = "estimate", aggregation = "mean", mcse = "sd"))
)

pos_prob_metric(name = "pos_prob")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "pos_prob".

## Value

A `PosProbMetric` object.

A `PosProbMetric` object.

## Examples

``` r
pos_prob_metric()
#> <bayesim::PosProbMetric>
#>  @ name        : chr "pos_prob"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ schema      :List of 2
#>  .. $ mean    :List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
#>  .. $ by_param:List of 3
#>  ..  ..$ role       : chr "estimate"
#>  ..  ..$ aggregation: chr "mean"
#>  ..  ..$ mcse       : chr "sd"
```
