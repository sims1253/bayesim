# Coverage Metric

Coverage of true parameter values within credible intervals.

Constructor for CoverageMetric.

## Usage

``` r
CoverageMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  prob = 0.95
)

coverage_metric(name = "coverage", prob = 0.95)
```

## Arguments

- name:

  Character string naming the metric. Defaults to "coverage".

- needs:

  Character vector of required capabilities. Defaults to empty character
  vector.

- required:

  Logical indicating if metric failure causes task failure. Defaults to
  FALSE.

- prob:

  Numeric probability for the credible interval width. Defaults to 0.95.

## Value

A CoverageMetric object.
