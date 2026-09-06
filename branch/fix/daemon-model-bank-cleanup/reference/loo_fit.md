# Compute LOO-CV

Compute leave-one-out cross-validation using Pareto-smoothed importance
sampling (PSIS-LOO). Named `loo_fit` to avoid clashing with
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

## Usage

``` r
loo_fit(fitter, fit_result, log_lik = NULL)
```

## Arguments

- fitter:

  An S7 Fitter object

- fit_result:

  A `bayesim_fit_result` object from
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)

- log_lik:

  Optional pointwise log-likelihood matrix (S x N, draws x observations)
  for the training set, as returned by
  [`log_lik_matrix()`](https://sims1253.github.io/bayesim/reference/log_lik_matrix.md).
  When supplied, methods should use it instead of recomputing their own;
  [`build_loo_context()`](https://sims1253.github.io/bayesim/reference/build_loo_context.md)
  passes the matrix it already computed so the weighted-prediction
  (PSIS) path pays for it once per task (#73). NULL (the default, and
  for standalone calls) means the method computes its own.

## Value

A list containing:

- `elpd`: Expected log predictive density (scalar)

- `p_loo`: Effective number of parameters (scalar)

- `elpd_se`: Standard error of ELPD (scalar)

- `pareto_k`: Pareto k diagnostic values (vector of length N)

- `r_eff`: Chain-aware relative efficiencies used for the summary
  (vector of length N), or NULL when none were computed (e.g. i.i.d.
  draws).
  [`build_loo_context()`](https://sims1253.github.io/bayesim/reference/build_loo_context.md)
  reuses it for the PSIS object (#73).

- Additional loo-specific diagnostics
