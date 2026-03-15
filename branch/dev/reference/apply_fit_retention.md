# Apply retention policy to fit result

Removes fields from fit result based on retention policy to reduce
memory footprint.

## Usage

``` r
apply_fit_retention(fit_result, retain, data_bundle = NULL)
```

## Arguments

- fit_result:

  A bayesim_fit_result object

- retain:

  Character vector of retention options specifying what to keep

- data_bundle:

  Optional data bundle; if provided and "data" not in retain, removes
  data_bundle from fit_result

## Value

Modified bayesim_fit_result object with non-retained fields removed
