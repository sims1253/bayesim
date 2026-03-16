# Posterior Mean Metric

Posterior mean estimates for parameters.

Constructor for PosteriorMeanMetric.

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
