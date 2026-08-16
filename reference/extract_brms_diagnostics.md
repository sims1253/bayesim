# Extract brms diagnostics

Computes rhat/ESS extrema over **all** parameters (fixed, group-level,
distributional, sigma) via
[`posterior::summarise_draws`](https://mc-stan.org/posterior/reference/draws_summary.html),
not just the fixed effects from `summary(fit)` (A3). Divergences and
max-treedepth hits come from the sampler diagnostics. `lp__` is excluded
as it is not a parameter of interest.

## Usage

``` r
extract_brms_diagnostics(fit)
```

## Arguments

- fit:

  brms fit object

## Value

Named list of diagnostics with `rhat_max`, `ess_bulk_min`,
`ess_tail_min`, `divergent`, `max_treedepth`.
