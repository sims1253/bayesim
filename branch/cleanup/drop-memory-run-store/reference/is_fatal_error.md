# Check if an error is fatal

Checks whether an error condition is fatal (should stop the entire
simulation run). Fatal errors include configuration errors, contract
violations, checkpoint failures, and internal consistency errors.

## Usage

``` r
is_fatal_error(cond)
```

## Arguments

- cond:

  A condition object to test

## Value

TRUE if the error is fatal, FALSE otherwise

## Examples

``` r
if (FALSE) { # \dontrun{
tryCatch(
  bayesim_config_error("Invalid parameter"),
  bayesim_error = function(cond) {
    if (is_fatal_error(cond)) {
      # Stop everything
    }
  }
)
} # }
```
