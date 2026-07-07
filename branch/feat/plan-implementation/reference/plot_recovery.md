# Plot parameter recovery (truth vs posterior estimate)

Scatter of posterior-mean estimates against true parameter values, per
task, with credible-interval segments. Faceted by a condition column
when `by` is supplied. Requires `posterior_summary_metric` to have been
computed and the truth recorded (E1).

## Usage

``` r
plot_recovery(result, var, by = NULL)
```

## Arguments

- result:

  A `bayesim_simulation_result`.

- var:

  Parameter name (a vars_of_interest entry).

- by:

  Optional name of a condition column to facet by (E7).

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_recovery(result, "b_x")
plot_recovery(result, "b_x", by = "data_n")
} # }
```
