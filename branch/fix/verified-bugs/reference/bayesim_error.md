# Base Bayesim Error

Base class for all bayesim errors. All bayesim-specific errors inherit
from this class in addition to "error" and "condition".

## Usage

``` r
bayesim_error(message, call = NULL)
```

## Arguments

- message:

  The error message

- call:

  The call that caused the error (optional)

## Value

An error condition object
