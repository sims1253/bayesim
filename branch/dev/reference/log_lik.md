# Compute Pointwise Log-Likelihood

Compute pointwise log-likelihood values.

## Usage

``` r
log_lik(fitter, fit_result, newdata = NULL)
```

## Arguments

- fitter:

  An S7 Fitter object

- fit_result:

  A `bayesim_fit_result` object from
  [`fit()`](https://sims1253.github.io/bayesim/reference/fit.md)

- newdata:

  Data frame with observations. If NULL, uses training data.

## Value

A matrix with dimensions N x S where:

- N = number of observations

- S = number of posterior draws

- Entry (i, s) is log p(y_i \| parameters_s)
