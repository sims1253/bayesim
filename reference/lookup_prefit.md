# Look up a prefit in the model bank

Computes the spec hash and returns the matching prefit object, or NULL
if the spec is not in the bank (caller falls back to a fresh compile).

## Usage

``` r
lookup_prefit(model_bank, formula, family, prior, stanvars, backend)
```

## Arguments

- model_bank:

  A named list of prefit objects (from
  [`build_model_bank()`](https://sims1253.github.io/bayesim/reference/build_model_bank.md)),
  or NULL.

- formula:

  A formula, brmsformula, or string.

- family:

  A brms family object, string, or NULL.

- prior:

  A brms prior object or NULL.

- stanvars:

  A brms stanvars object or NULL.

- backend:

  Character scalar.

## Value

A `brmsfit` prefit object, or NULL.
