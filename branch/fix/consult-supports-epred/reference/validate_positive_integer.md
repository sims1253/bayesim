# Validate a Positive Integer

Returns an S7 property validator function that requires `value` to be a
single non-NA integer greater than or equal to 1. The returned function
returns `NULL` when valid, otherwise the string `"<message>"` (defaults
to `"must be a positive integer"`).

## Usage

``` r
validate_positive_integer(message = "must be a positive integer")
```

## Arguments

- message:

  Character scalar; the error string returned when `value` is invalid.
  Defaults to `"must be a positive integer"`. Pass a
  property-name-prefixed string to preserve a previously established
  error message verbatim.

## Value

A `function(value)` suitable as an S7 property `validator`.
