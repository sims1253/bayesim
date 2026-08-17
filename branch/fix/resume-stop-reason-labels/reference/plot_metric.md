# Plot a metric across conditions

General-purpose plot of a metric column (or its aggregated mean) against
a condition variable, with optional facets. Requires ggplot2.

## Usage

``` r
plot_metric(result, metric, x = NULL, facets = NULL)
```

## Arguments

- result:

  A `bayesim_simulation_result` or summary tibble.

- metric:

  Name of the metric column to plot.

- x:

  Optional x-axis variable (condition). Defaults to `task_id`.

- facets:

  Optional facet variable.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_metric(result, "rmse__value", x = "n")
} # }
```
