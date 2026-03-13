# Validate Metric Interface

Validates that a metric object conforms to the Metric S7 class
interface. The metric must be an S7 object that inherits from the Metric
class and implements the compute method.

## Usage

``` r
validate_metric_interface(metric)
```

## Arguments

- metric:

  An S7 object to validate as a Metric.

## Value

The input `metric`, invisibly, if validation passes.

## Details

The metric must satisfy the following requirements:

- Must be an S7 object (checked via S7::S7_inherits())

- Must inherit from the "Metric" class

- Must implement the
  `compute(fit_result, data_bundle, context, task_ctx)` method

- Must have a `name` property that is a non-empty character string

## Errors

Throws a `bayesim_contract_error` condition if validation fails.

## See also

[Metric](https://sims1253.github.io/bayesim/reference/Metric.md),
[`validate_metric_output()`](https://sims1253.github.io/bayesim/reference/validate_metric_output.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a custom metric using proper S7 syntax
MyMetric <- S7::new_class(
  "MyMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "my_metric")
  )
)
# Register the compute method separately
S7::method(compute, MyMetric) <- function(metric, fit_result, data_bundle, context, task_ctx) {
  list(value = 1.0)
}
my_metric <- MyMetric(name = "my_metric")
validate_metric_interface(my_metric)
} # }
```
