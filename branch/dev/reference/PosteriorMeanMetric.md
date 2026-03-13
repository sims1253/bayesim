# Posterior Mean Metric

Posterior mean estimates for parameters.

`posterior_mean_metric()` is a constructor function that creates a
PosteriorMeanMetric instance with appropriate defaults. Use this
constructor to work around S7's property default inheritance issue.

## Usage

``` r
PosteriorMeanMetric(
  name = character(0),
  needs = character(0),
  required = FALSE
)

posterior_mean_metric(name = "posterior_mean")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "posterior_mean".

- needs:

  Character vector of required capabilities. Defaults to empty character
  vector.

- required:

  Logical indicating if metric failure causes task failure. Defaults to
  FALSE.

## Value

A PosteriorMeanMetric object.

## Examples

``` r
posterior_mean_metric()
#> <bayesim::PosteriorMeanMetric>
#>  @ name    : chr "posterior_mean"
#>  @ needs   : chr(0) 
#>  @ required: logi FALSE
posterior_mean_metric(name = "my_posterior_mean")
#> <bayesim::PosteriorMeanMetric>
#>  @ name    : chr "my_posterior_mean"
#>  @ needs   : chr(0) 
#>  @ required: logi FALSE
```
