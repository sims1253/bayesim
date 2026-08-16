# Checkpoint read/write error (Fatal)

Raised when checkpoint file operations fail, such as reading corrupted
checkpoint data or failing to write checkpoint files. This is a fatal
error that stops the entire simulation run.

## Usage

``` r
bayesim_checkpoint_error(message, call = NULL)
```

## Arguments

- message:

  The error message

- call:

  The call that caused the error (optional)

## Value

An error condition object with class c("bayesim_checkpoint_error",
"bayesim_error", "error", "condition")

## See also

[`is_fatal_error()`](https://sims1253.github.io/bayesim/reference/is_fatal_error.md)
