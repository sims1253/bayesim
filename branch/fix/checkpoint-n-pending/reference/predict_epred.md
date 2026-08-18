# Compute Posterior Expectation Predictions

Compute posterior draws of the expected value of the response
distribution (mu, without observation noise). This is the `epred`
quantity used by brms' `loo_R2()` and bayesim's
[`r2_loo_metric()`](https://sims1253.github.io/bayesim/reference/R2LooMetric.md)
(F3). It must NOT include observation-level noise — only the model's
conditional mean.

Fitters that cannot provide expectation predictions should return
`NULL`; the `r2_loo` metric then degrades to NA. The default Fitter
method returns NULL.

## Usage

``` r
predict_epred(fitter, fit_result, newdata = NULL)
```

## Arguments

- fitter:

  An S7 Fitter object

- fit_result:

  A `bayesim_fit_result` object from
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)

- newdata:

  Data frame with observations. If NULL, uses training data.

## Value

A matrix with dimensions S x N (draws x observations), or NULL if not
supported by the fitter.
