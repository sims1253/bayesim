# Run a fresh brms fit (no model bank)

Shared helper for the fallback path (`precompile = FALSE` or bank miss).

## Usage

``` r
run_fresh_brms(
  fitter,
  data_bundle,
  fit_spec,
  seed,
  formula,
  family,
  prior,
  stanvars
)
```
