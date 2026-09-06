# Validate a Metric Object

Validates that a metric object is an S7 instance of the Metric class
with a valid `name` property. This is the canonical metric validator (B3
merges the former internal `validate_metric_interface` into this
exported name). Method existence is not checked because S7 dispatches
via generics and the base class raises errors for unimplemented abstract
methods.

When representative values are supplied (`fit_result` plus
`data_bundle`), validation additionally executes
[`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md)
once and checks that its output passes the metric output schema and that
every field declared in the metric's `schema` is produced. This catches
broken external metrics before a run starts instead of as thousands of
per-task failures.

## Usage

``` r
validate_metric(
  metric,
  fit_result = NULL,
  data_bundle = NULL,
  context = NULL,
  task_ctx = NULL
)
```

## Arguments

- metric:

  An S7 object to validate as a Metric.

- fit_result:

  Optional representative `bayesim_fit_result`
  ([`new_fit_result()`](https://sims1253.github.io/bayesim/reference/new_fit_result.md))
  driving a conformance execution of
  [`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md).
  Must be supplied together with `data_bundle`.

- data_bundle:

  Optional representative data bundle for the conformance execution.

- context:

  Optional representative context list (e.g. predictions, log_lik).
  Defaults to an empty list; metrics should degrade to `NA` rather than
  erroring when a needed context element is missing.

- task_ctx:

  Optional representative task context. Defaults to a stable
  placeholder.

## Value

The input `metric`, invisibly, if validation passes.

## Errors

Throws a `bayesim_contract_error` condition if validation fails (or its
subclass `bayesim_validation_error` for representative-execution
failures).

## See also

[Metric](https://sims1253.github.io/bayesim/reference/Metric.md),
[`validate_metric_output()`](https://sims1253.github.io/bayesim/reference/validate_metric_output.md),
[`validate_fitter()`](https://sims1253.github.io/bayesim/reference/validate_fitter.md)
