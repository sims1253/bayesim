# Apply retention policy to task result

Adds optional retained fields from the fit result and data bundle to the
task result based on the retention policy.

## Usage

``` r
apply_task_retention(task_result, fit_result, data_bundle, retain)
```

## Arguments

- task_result:

  A bayesim_task_result object to modify

- fit_result:

  A bayesim_fit_result object (before retention applied)

- data_bundle:

  Data bundle from generator containing train/test data

- retain:

  Character vector of retention options specifying what to keep

## Value

Modified bayesim_task_result object with retained fields added
