# Resolve retention specification to context-specific form

Converts a retention profile name, character vector, or named list into
a canonical named list with success/warning/error entries.

## Usage

``` r
resolve_retention_spec(retain)
```

## Arguments

- retain:

  Character vector/profile or named list with retention options

## Value

Named list with entries for "success", "warning", and "error" contexts
