# Resolve requested cleaned var names to actual draws-matrix columns

Maps each name in `vars` to itself, then to `b_<name>`, in the set of
`draws_colnames`, mirroring the both-directions lookup
[`.extract_truth()`](https://sims1253.github.io/bayesim/reference/dot-extract_truth.md)
already performs. This fixes the F2 mismatch where generators strip the
`b_` prefix (so `vars_of_interest` holds cleaned names like
`c("x", "Intercept")`) but brms draws matrices keep it
(`c("b_x","b_Intercept", "sigma")`).

## Usage

``` r
resolve_draw_columns(vars, draws_colnames)
```

## Arguments

- vars:

  Character vector of cleaned names (e.g.
  `c("x","Intercept", "sigma")`).

- draws_colnames:

  Character vector of available draws-matrix column names.

## Value

A **named** character vector: names are the cleaned `vars` (use these
for output field naming), values are the actual draws column to read.
Empty input returns `character(0)` (names preserved).

## Details

Errors with a `bayesim_config_error` (same condition class as
[`.extract_truth()`](https://sims1253.github.io/bayesim/reference/dot-extract_truth.md))
when a requested var is genuinely absent. A silent NA would corrupt SBC
ranks and credible intervals without diagnostics.
