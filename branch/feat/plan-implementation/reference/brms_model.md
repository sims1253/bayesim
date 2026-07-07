# Construct a single brms model specification

Builds a one-row list of brms model components (formula, family, prior,
stanvars, stan_file) suitable for assembling a fit_grid via
[`model_grid()`](https://sims1253.github.io/bayesim/reference/model_grid.md).
Inputs are validated at construction so errors surface before
compilation.

## Usage

``` r
brms_model(
  formula,
  family = NULL,
  prior = NULL,
  stanvars = NULL,
  stan_file = NULL
)
```

## Arguments

- formula:

  A formula or brmsformula.

- family:

  A brms family (e.g.
  [`gaussian()`](https://rdrr.io/r/stats/family.html), `student()`), or
  NULL for the brms default.

- prior:

  A brms prior object, or NULL.

- stanvars:

  A brms `stanvars` object, or NULL.

- stan_file:

  Optional path to a `.stan` file (passed through to brms).

## Value

A named list with elements `formula`, `family`, `prior`, `stanvars`,
`stan_file`.

## See also

[`model_grid()`](https://sims1253.github.io/bayesim/reference/model_grid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
brms_model(y ~ x, family = gaussian())
} # }
```
