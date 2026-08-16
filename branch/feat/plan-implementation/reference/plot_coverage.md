# Plot coverage rates per condition/parameter

Point-range plot of credible-interval coverage with MCSE error bars (E7
redesign: was a bar plot of a continuous coverage_mean). Coverage and
its MCSE come from
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md);
each point is a condition x estimand cell, with a dashed reference line
at the nominal rate. Requires
[`posterior_summary_metric()`](https://sims1253.github.io/bayesim/reference/PosteriorSummaryMetric.md)
(and recorded truths, E1).

## Usage

``` r
plot_coverage(result, nominal = NULL, by = NULL)
```

## Arguments

- result:

  A `bayesim_simulation_result`.

- nominal:

  Nominal coverage rate. Defaults to the interval probability recorded
  by the metric schema, or 0.95 for legacy results.

- by:

  Character vector of condition columns (passed to
  [`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_coverage(result)
} # }
```
