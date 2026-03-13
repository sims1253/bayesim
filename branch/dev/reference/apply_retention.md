# Apply retention policy to fit result

Removes large objects from the fit result based on the retention policy
to manage memory usage during simulation runs.

## Usage

``` r
apply_retention(fit_result, data_bundle, retain)
```

## Arguments

- fit_result:

  A bayesim_fit_result object

- data_bundle:

  A data bundle list (currently unused, for future expansion)

- retain:

  Character vector specifying what to retain. Options:

  - "fit": Keep the raw fit object

  - "draws": Keep the posterior draws matrix

  - "diagnostics": Keep the diagnostics list

## Value

The fit_result with specified elements set to NULL based on the
retention policy.

## Details

By default, only metrics and diagnostics are retained. Large objects
like the raw fit and draws matrix are removed to minimize memory usage.
Users can override this by including "fit" or "draws" in the retain
vector.
