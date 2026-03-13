# Apply retention policy to fit result

Removes fields from fit result based on retention policy to reduce
memory footprint. This function modifies the fit_result in place by
setting unwanted fields to NULL.

## Usage

``` r
apply_fit_retention(fit_result, retain)
```

## Arguments

- fit_result:

  A bayesim_fit_result object

- retain:

  Character vector of retention options specifying what to keep

## Value

Modified bayesim_fit_result object with non-retained fields removed
