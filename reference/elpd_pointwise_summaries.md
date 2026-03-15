# Pointwise elps summaries

Convenience function to collect quantiles and summaries of the pointwise
elpd estimates instead of just the main estimates. The returned
summaries are all sample size independent.

## Usage

``` r
elpd_pointwise_summaries(fit, quantiles, newdata = NULL)
```

## Arguments

- fit:

  A brmsfit object.

- quantiles:

  A vector of quantiles of interest.

- newdata:

  If supplied, returns the summaries for
  [`elpd_test()`](https://sims1253.github.io/bayesim/reference/elpd_test.md)
  otherwise, returns `elpd` summaries by default.

## Value

A named list of summaries.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- brms::brm(y ~ 1, data = list(rnorm(1000)))
elpd_pointwise_summaries(fit, seq(0.1, 0.9, length.out = 9))
elpd_pointwise_summaries(fit, seq(0.1, 0.9, length.out = 9), rnorm(1000))
} # }
```
