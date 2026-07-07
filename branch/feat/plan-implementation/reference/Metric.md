# Metric Abstract Class

Abstract base class for defining simulation metrics in bayesim. Metrics
compute summary statistics or diagnostic values from fitted model
results and data bundles.

## Usage

``` r
Metric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean"
)
```

## Arguments

- name:

  Character identifier for the metric. Used as a prefix when flattening
  metric output to column names.

- needs:

  Character vector of required capabilities from the fitter. Common
  values include "predictions", "log_lik", "loo". The metric will only
  receive these values in the context if the fitter provides them.

- required:

  Logical indicating whether metric failure causes task failure. If
  TRUE, an error in computing this metric will propagate and fail the
  entire task. If FALSE (default), metric failure results in NA values
  being recorded.

- summary_type:

  Character; how
  [`summarize_simulation()`](https://sims1253.github.io/bayesim/reference/summarize_simulation.md)
  aggregates this metric's flattened columns: `"mean"` (default,
  sd/sqrt(n) MCSE), `"proportion"` (coverage-style sqrt(p(1-p)/n) MCSE),
  or `"none"`.

## Value

An S7 class object representing a Metric.

## Methods

The `compute()` S7 generic must be implemented by subclasses.

- compute(metric, fit_result, data_bundle, context, task_ctx):

  Compute metric values from a fitted model result. This method must be
  implemented by subclasses.

      \itemize{
        \item metric: The Metric S7 object
        \item fit_result: A bayesim_fit_result object containing the fitted
          model output (draws, diagnostics, etc.)
        \item data_bundle: A list containing data-related objects including
          train (training data), test (test data if applicable), response
          (response variable), true_params (true parameter values if known),
        \item context: A list with precomputed values based on the metric's
          `needs` property. May include predictions, log_lik (log-likelihood
          values), loo (leave-one-out cross-validation results), etc.
        \item task_ctx: A list with task identification information including
          task_id, data_idx, fit_idx, rep_idx for tracking and debugging.
      }

      Returns: A named list with metric values. Names must be non-empty strings.
      Values must be one of:
      \itemize{
        \item scalar atomic (logical, integer, double, character)
        \item named numeric vector
      }

      No nested data frames or matrices are allowed in the output.

## Metric Output Schema

The compute method must return a named list conforming to the following
schema:

- All elements must have non-empty names

- Values must be scalar atomic types or named numeric vectors

- No nested data frames or matrices allowed

- The engine flattens output with prefix `<metric_name>__<field>`

## Examples

``` r
# Define a custom RMSE metric
RMSEMetric <- S7::new_class(
  "RMSEMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "rmse"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

S7::method(compute_metric, RMSEMetric) <- function(
  metric, fit_result, data_bundle, context, task_ctx
) {
  preds <- context$predictions
  actual <- data_bundle$test[[data_bundle$response]]
  list(
    value = sqrt(mean((preds - actual)^2)),
    n_obs = length(actual)
  )
}
```
