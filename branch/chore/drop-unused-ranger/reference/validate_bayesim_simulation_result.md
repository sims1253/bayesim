# Validate a bayesim_simulation_result object

Performs consistency checks on a bayesim_simulation_result object.

## Usage

``` r
validate_bayesim_simulation_result(x)
```

## Arguments

- x:

  A bayesim_simulation_result object to validate

## Value

The input object, invisibly, if validation passes

## Errors

Throws an error if validation fails, with a message indicating the
specific validation problem (e.g., class mismatch, invalid
config_fingerprint, invalid task_results elements, non-data.frame
summary/errors, etc.).
