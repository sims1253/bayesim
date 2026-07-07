# Per-observation relative efficiency from the fit's chain structure

Mirrors brms' internal `r_eff_log_lik()` by deriving the chain id per
draw from `posterior::as_draws_df(fit)$.chain` and passing it to
`loo::relative_eff(exp(ll), chain_id = ...)`. `ll` here is S x N (draws
x observations, matching brms' log_lik orientation): `relative_eff`
wants the draws dimension as rows so chain_id (length S) matches
`nrow(ll)`. Returns one r_eff per observation (length N). Returns NULL
when chain information is unavailable (e.g. MockFitter or a fitter whose
fit lacks a `.chain` variable), in which case the caller falls back to
`r_eff = NULL`.

## Usage

``` r
relative_eff_from_chains(fitter, fit_result, ll)
```
