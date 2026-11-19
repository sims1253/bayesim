# Compute Metric Values

Compute metric values from a fitted model result. This generic must be
implemented by Metric subclasses. Named `compute_metric` (rather than
`compute`) so that bayesim does not mask
[dplyr::compute](https://dplyr.tidyverse.org/reference/compute.html).

## Usage

``` r
compute_metric(metric, fit_result, data_bundle, context, task_ctx)
```

## Arguments

- metric:

  An S7 Metric object

- fit_result:

  A bayesim_fit_result object containing the fitted model output

- data_bundle:

  A list containing data-related objects

- context:

  A list with precomputed values based on the metric's `needs`

- task_ctx:

  A list with task identification information

## Value

A named list with metric values conforming to the metric output schema
