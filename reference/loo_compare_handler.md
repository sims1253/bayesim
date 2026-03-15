# Compare LOO objects

Compare LOO objects

## Usage

``` r
loo_compare_handler(loo_object_matrix, predictive_metrics)
```

## Arguments

- loo_object_matrix:

  List of LOO objects

- predictive_metrics:

  Character vector of metrics to compare

## Value

A data frame with delta and se_delta for each metric

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package and fitted models
loo_compare_handler(loo_objects, metrics)
} # }
```
