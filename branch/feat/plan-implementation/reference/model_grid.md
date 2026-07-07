# Assemble a tibble of model specifications for a fit_grid

Combines named
[`brms_model()`](https://sims1253.github.io/bayesim/reference/brms_model.md)
specs into a single data frame with one row per model and
`formula`/`family`/`prior`/`stanvars`/`stan_file` list-columns, ready to
pass as `fit_grid` to
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md).
A `model` column holds the spec names and lands in the summary as
`fit_model`.

## Usage

``` r
model_grid(...)
```

## Arguments

- ...:

  Named
  [`brms_model()`](https://sims1253.github.io/bayesim/reference/brms_model.md)
  specs (or any named lists with the same component names).

## Value

A tibble with columns `model`, `formula`, `family`, `prior`, `stanvars`,
`stan_file`.

## See also

[`brms_model()`](https://sims1253.github.io/bayesim/reference/brms_model.md),
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)

## Examples

``` r
if (FALSE) { # \dontrun{
grid <- model_grid(
  gaussian = brms_model(y ~ x, gaussian()),
  student  = brms_model(y ~ x, student())
)
} # }
```
