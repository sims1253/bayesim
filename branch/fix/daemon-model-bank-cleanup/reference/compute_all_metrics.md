# Compute all metrics for a task

Iterates through all metrics and computes their values, handling errors
gracefully based on each metric's `required` property.

## Usage

``` r
compute_all_metrics(
  fit_result,
  data_bundle,
  context,
  metrics,
  task_ctx,
  result_path = NULL
)
```

## Arguments

- fit_result:

  A bayesim_fit_result object from a successful fit

- data_bundle:

  A validated data bundle list

- context:

  A list with precomputed values (from build_metric_context)

- metrics:

  List of S7 Metric objects

- task_ctx:

  Task context with identification information

## Value

A named list with:

- `metrics`: Named list of computed and flattened metric values

- `warnings`: Character vector of any warning messages

## Details

For each metric:

- If the metric has `required = TRUE` and fails, the error is re-thrown

- If the metric has `required = FALSE` and fails, NA values are returned
  with error information stored in `<metric_name>__error_class` and
  `<metric_name>__error_message`

- All metric outputs are flattened using
  [`flatten_metric_output()`](https://sims1253.github.io/bayesim/reference/flatten_metric_output.md)
