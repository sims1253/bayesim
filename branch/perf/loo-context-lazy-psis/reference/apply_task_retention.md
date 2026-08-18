# Apply retention policy to task result

Adds optional retained fields from the fit result and data bundle to the
task result based on the retention policy.

## Usage

``` r
apply_task_retention(
  task_result,
  fit_result,
  data_bundle,
  retain,
  predictions = NULL
)
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

- predictions:

  Prediction context returned by
  [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md),
  or NULL.

## Value

Modified bayesim_task_result object with retained fields added
