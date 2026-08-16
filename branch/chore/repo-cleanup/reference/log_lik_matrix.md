# Compute Pointwise Log-Likelihood

Compute pointwise log-likelihood values. Named `log_lik_matrix` (rather
than `log_lik`) so that bayesim does not mask
[brms::log_lik](https://mc-stan.org/rstantools/reference/log_lik.html)
or the rstantools `log_lik` generic for users who load bayesim alongside
brms.

## Usage

``` r
log_lik_matrix(fitter, fit_result, newdata = NULL)
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

A matrix with dimensions S x N where:

- S = number of posterior draws (rows)

- N = number of observations (columns)

- Entry (s, i) is log p(y_i \| parameters_s)

This is the brms/loo convention (draws as rows), matching
[`brms::log_lik`](https://mc-stan.org/rstantools/reference/log_lik.html)
and what
[`loo::psis`](https://mc-stan.org/loo/reference/psis.html)/[`loo::E_loo`](https://mc-stan.org/loo/reference/E_loo.html)/[`loo::relative_eff`](https://mc-stan.org/loo/reference/relative_eff.html)
expect.
