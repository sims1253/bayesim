# Compute LOO-CV

Compute leave-one-out cross-validation using Pareto-smoothed importance
sampling (PSIS-LOO). Named `loo_fit` to avoid clashing with
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

## Usage

``` r
loo_fit(fitter, fit_result)
```

## Arguments

- fitter:

  An S7 Fitter object

- fit_result:

  A `bayesim_fit_result` object from
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)

## Value

A list containing:

- `elpd`: Expected log predictive density (scalar)

- `p_loo`: Effective number of parameters (scalar)

- `elpd_se`: Standard error of ELPD (scalar)

- `pareto_k`: Pareto k diagnostic values (vector of length N)

- Additional loo-specific diagnostics
