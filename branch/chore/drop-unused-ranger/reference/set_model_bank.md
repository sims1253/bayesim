# Set the session model bank

Stores the model bank in `options(bayesim.model_bank)`. Called on the
controller after the bank is built, and pushed to each daemon via
[`mirai::everywhere()`](https://mirai.r-lib.org/reference/everywhere.html).
Pass NULL to clear it.

## Usage

``` r
set_model_bank(bank)
```

## Arguments

- bank:

  A named list of `brmsfit` prefit objects keyed by
  [`model_spec_hash()`](https://sims1253.github.io/bayesim/reference/model_spec_hash.md),
  or NULL.

## Value

Invisible NULL. Called for its side effect.
