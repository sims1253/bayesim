# Build the LOO context (elpd summary + PSIS object + epred)

Constructs the full LOO context for F3's rmse_loo / r2_loo metrics.
Computes the elpd/p_loo/pareto_k summary (as the legacy
[`loo_fit()`](https://sims1253.github.io/bayesim/reference/loo_fit.md)
did), the PSIS object (for
[`loo::E_loo()`](https://mc-stan.org/loo/reference/E_loo.html) weighted
predictions), the pointwise log-likelihood matrix, and the posterior
expectation predictions (epred) — all once, shared across metrics.

## Usage

``` r
build_loo_context(fitter, fit_result)
```

## Value

A list with elements `loo`, `psis`, `log_lik`, `epred`, or NULL on
failure. `psis`/`log_lik`/`epred` may be individually NULL if
unavailable.

## Details

The PSIS object uses `loo::psis(-ll, r_eff)` with per-observation
relative-efficiency factors derived from the fit's chain structure via
`posterior::as_draws_df(fit)$.chain` (matches brms' internal
`r_eff_log_lik` exactly). Falls back to `r_eff = NULL` (with a captured
warning) when chain structure is unavailable, which is mathematically
valid but slightly less accurate.

epred must be the posterior expectation (mu, no observation noise); for
brms this is
[`brms::posterior_epred`](https://mc-stan.org/rstantools/reference/posterior_epred.html).
Only fitters with `supports_epred = TRUE` are asked for it via
[`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md);
otherwise epred is NULL and the consuming metrics (r2_loo, rmse_loo)
degrade to NA.
