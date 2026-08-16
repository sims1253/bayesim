# Build context for metrics computation

Precomputes predictions, log_lik, and loo based on metric needs and
retained predictions. Only computes requested context that the fitter
supports.

## Usage

``` r
build_metric_context(
  fit_result,
  fitter,
  data_bundle,
  metrics,
  seed = NULL,
  retain = character()
)
```

## Arguments

- fit_result:

  A bayesim_fit_result object from a successful fit

- fitter:

  S7 Fitter object

- data_bundle:

  A validated data bundle list

- metrics:

  List of S7 Metric objects

- retain:

  Retention specification. Requesting `"predictions"` computes
  predictions even when no metric needs them.

## Value

A named list containing any of:

- `predictions`: Prediction results from the fitter

- `log_lik`: Pointwise log-likelihood matrix (S x N)

- `loo`: LOO-CV results

## Details

The function inspects the `needs` property of each metric to determine
what context elements are required. It then checks if the fitter
supports each capability and computes them. Any errors during
computation result in NULL values for that context element.

Evaluation data: `predictions` and `log_lik` are computed on the TEST
set when `data_bundle$test` is present, otherwise on the training set.
Every built-in metric that consumes them (`pred_*`, `elpd_test`,
`r2_test`) compares against the test response, so the predictions must
be for the test rows. The LOO context is always built on the training
set — leave-one-out is in-sample by construction.
