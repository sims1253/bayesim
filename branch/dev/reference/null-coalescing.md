# Null coalescing operator

Returns `x` if not NULL, otherwise `y`.

## Usage

``` r
x %||% y
```

## Arguments

- x:

  Value to check

- y:

  Default value if x is NULL

## Value

`x` if not NULL, otherwise `y`

## Examples

``` r
NULL %||% "default"  # returns "default"
#> [1] "default"
"value" %||% "default"  # returns "value"
#> [1] "value"
```
