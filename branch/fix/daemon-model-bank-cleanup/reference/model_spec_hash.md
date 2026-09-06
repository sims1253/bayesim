# Compute a stable hash for a model spec

Produces a stable string key identifying a brms model spec up to the
fields that affect Stan code generation: formula, family, prior,
stanvars, and backend. Uses
[`digest::digest()`](https://eddelbuettel.github.io/digest/man/digest.html)
on a normalized list representation, which is deterministic across R
sessions (unlike environment-bound object hashing).

## Usage

``` r
model_spec_hash(formula, family, prior, stanvars, backend)
```

## Arguments

- formula:

  A formula, brmsformula, or string (resolved to brmsformula).

- family:

  A brms family object, string, or NULL (resolved).

- prior:

  A brms prior object or NULL.

- stanvars:

  A brms stanvars object or NULL.

- backend:

  Character scalar, e.g. "cmdstanr".

## Value

A length-1 character string hash.
