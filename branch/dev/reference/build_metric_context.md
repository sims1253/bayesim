# Build context for metrics computation

Precomputes predictions, log_lik, loo based on metric needs. Only
computes what is needed by the metrics, and only if the fitter supports
it.

## Usage

``` r
build_metric_context(fit_result, fitter, data_bundle, metrics, seed = NULL)
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

## Value

A named list containing any of:

- `predictions`: Prediction results from the fitter

- `log_lik`: Pointwise log-likelihood matrix

- `loo`: LOO-CV results

## Details

The function inspects the `needs` property of each metric to determine
what context elements are required. It then checks if the fitter
supports each capability and computes them. Any errors during
computation result in NULL values for that context element.
