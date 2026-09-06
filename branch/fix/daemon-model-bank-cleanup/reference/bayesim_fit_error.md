# Model fitting error (Recoverable)

Raised when model fitting fails for a specific task, such as convergence
failures or numerical issues. This is a task-level recoverable error
that allows the simulation to continue with other tasks.

## Usage

``` r
bayesim_fit_error(message, call = NULL)
```

## Arguments

- message:

  The error message

- call:

  The call that caused the error (optional)

## Value

An error condition object with class c("bayesim_fit_error",
"bayesim_error", "error", "condition")

## See also

[`is_recoverable_error()`](https://sims1253.github.io/bayesim/reference/is_recoverable_error.md)
