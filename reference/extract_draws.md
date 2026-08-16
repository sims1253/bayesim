# Extract Posterior Draws

Extract posterior draws from a fitted model.

## Usage

``` r
extract_draws(fitter, fit_result, variables = NULL)
```

## Arguments

- fitter:

  An S7 Fitter object

- fit_result:

  A `bayesim_fit_result` object from
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)

- variables:

  Character vector of variable names to extract. If NULL, extracts all
  variables.

## Value

A matrix with dimensions S x P where:

- S = number of posterior draws (iterations x chains)

- P = number of parameters/variables

- Column names match variable names
