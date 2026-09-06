# Check if an error is recoverable (task-level)

Checks whether an error condition is recoverable at the task level.
Recoverable errors include data generation errors, model fitting errors,
and metric computation errors. The simulation can continue with other
tasks when these errors occur.

## Usage

``` r
is_recoverable_error(cond)
```

## Arguments

- cond:

  A condition object to test

## Value

TRUE if the error is recoverable, FALSE otherwise

## Examples

``` r
if (FALSE) { # \dontrun{
tryCatch(
  bayesim_fit_error("Model did not converge"),
  bayesim_error = function(cond) {
    if (is_recoverable_error(cond)) {
      # Log error and continue with next task
    }
  }
)
} # }
```
