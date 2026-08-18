# Flatten a nested list with prefix

Flattens a named list with nested named numeric vectors, using a
prefix\_\_name\_\_subname naming convention for the flattened elements.

## Usage

``` r
flatten_with_prefix(x, prefix)
```

## Arguments

- x:

  A named list, possibly containing nested named numeric vectors

- prefix:

  Character string to use as prefix for flattened names

## Value

A flattened named list where nested named numeric vectors are expanded.
When `prefix` is non-empty, flattened names use
"prefix\_\_name\_\_sub_name"; when `prefix` is empty (""), they use
"name\_\_sub_name". Scalar and unnamed elements are passed through
unchanged.

## Details

This function handles lists where some elements may be named numeric
vectors. Those vectors are flattened into individual scalar elements
with a double-underscore naming convention. The empty-prefix variant is
used internally by checkpointing code where the metric name already
serves as the outer namespace.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- list(a = 1, b = c(x = 2, y = 3), c = 4)
flatten_with_prefix(x, "param")
# Returns: list(a = 1, param__b__x = 2, param__b__y = 3, c = 4)
} # }
```
