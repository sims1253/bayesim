# Resolve retention specification

Converts retention profile name or explicit vector to canonical form. If
a profile name is provided, returns the corresponding options. If
explicit options are provided, validates them and returns only valid
ones.

## Usage

``` r
resolve_retention(retain)
```

## Arguments

- retain:

  Character vector of retention options or a profile name (one of
  "minimal", "standard", or "debug")

## Value

Character vector of valid retention options

## Examples

``` r
if (FALSE) { # \dontrun{
resolve_retention("minimal")
resolve_retention(c("metrics", "draws"))
} # }
```
