# Capture error information

Captures detailed information about an error condition, including class,
message, call, and a trimmed traceback.

## Usage

``` r
capture_error_info(e)
```

## Arguments

- e:

  A condition object captured by a stack-preserving handler or with an
  attached `rlang` trace.

## Value

A named list with elements:

- `error_class` - The class of the error

- `error_message` - The error message

- `call` - The call that caused the error (if available)

- `traceback` - Trimmed traceback (limited to 20 frames)
