# Extract real warmup/sample/total timings from a brmsfit

brms stores fit output as a `stanfit` object regardless of backend
(cmdstanr output is converted), so
[`rstan::get_elapsed_time()`](https://mc-stan.org/rstan/reference/stanfit-class.html)
works for both backends. Returns per-chain CPU seconds summed across
chains. Falls back to `NA_real_` for warmup/sample (and the provided
fallback total) when the stanfit timing is unavailable (e.g. a chains =
0 prefit, or an interrupted run).

## Usage

``` r
extract_brms_timings(fit, fallback_total)
```

## Arguments

- fit:

  A brmsfit object.

- fallback_total:

  Numeric length-1; used for `total` when timing can't be extracted
  (e.g. the wall-clock elapsed time from a timer).

## Value

Named list: `list(total, warmup, sample)`, all numeric scalars.
