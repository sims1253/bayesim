# Validate a bayesim_task_result object

Performs consistency checks on a bayesim_task_result object.

## Usage

``` r
validate_bayesim_task_result(x)
```

## Arguments

- x:

  A bayesim_task_result object to validate

## Value

The input object, invisibly, if validation passes

## Errors

Throws an error if validation fails, with a message indicating the
specific validation problem (e.g., class mismatch, invalid task_id or
status, missing metrics for successful tasks, missing error for failed
tasks, etc.).
