# Validate a bayesim_fit_result object

Performs consistency checks on a bayesim_fit_result object.

## Usage

``` r
validate_bayesim_fit_result(x)
```

## Arguments

- x:

  A bayesim_fit_result object to validate

## Value

The input object, invisibly, if validation passes

## Errors

Throws an error if validation fails, with a message indicating the
specific validation problem (e.g., class mismatch, success/error
inconsistency, invalid timing, missing draws colnames, etc.).
