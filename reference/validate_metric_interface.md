# Validate Metric Class and Name Property

Validates that a metric object is an S7 instance of the Metric class
with a valid `name` property. This is a lightweight check used
internally by
[`validate_simulation_config()`](https://sims1253.github.io/bayesim/reference/validate_simulation_config.md).
Method existence is not checked because S7 dispatches via generics and
the base class raises errors for unimplemented abstract methods.

## Usage

``` r
validate_metric_interface(metric)
```

## Arguments

- metric:

  An S7 object to validate as a Metric.

## Value

The input `metric`, invisibly, if validation passes.

## Errors

Throws a `bayesim_contract_error` condition if validation fails.

## See also

[Metric](https://sims1253.github.io/bayesim/reference/Metric.md),
[`validate_metric_output()`](https://sims1253.github.io/bayesim/reference/validate_metric_output.md)
