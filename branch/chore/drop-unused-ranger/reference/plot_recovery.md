# Plot parameter recovery (truth vs posterior estimate)

Scatter of posterior-mean estimates against true parameter values, per
task, with credible-interval segments. Faceted by a condition column
when `by` is supplied. Requires `posterior_summary_metric` to have been
computed and the truth recorded (E1).

## Usage

``` r
plot_recovery(result, estimand = NULL, by = NULL, var = NULL)
```

## Arguments

- result:

  A `bayesim_simulation_result`.

- estimand:

  Parameter name (a `vars_of_interest` entry); the preferred
  terminology, matching
  [`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md).

- by:

  Optional name of a condition column to facet by (E7).

- var:

  Legacy alias for `estimand`, kept for backward compatibility with
  earlier bayesim versions. It is a silent compatibility alias, not a
  deprecated argument, and emits no deprecation warning. When both
  `estimand` and `var` are supplied they must name the same parameter;
  conflicting values are an error.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_recovery(result, estimand = "b_x")
plot_recovery(result, estimand = "b_x", by = "data_n")
# Legacy spelling, identical result:
plot_recovery(result, var = "b_x")
} # }
```
