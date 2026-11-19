# Validation error (Fatal)

Raised when input validation fails. This is a fatal error that stops the
entire simulation run.

## Usage

``` r
bayesim_validation_error(message, call = NULL)
```

## Arguments

- message:

  The error message

- call:

  The call that caused the error (optional)

## Value

An error condition object with class c("bayesim_validation_error",
"bayesim_contract_error", "bayesim_error", "error", "condition")

## See also

[`is_fatal_error()`](https://sims1253.github.io/bayesim/reference/is_fatal_error.md)
