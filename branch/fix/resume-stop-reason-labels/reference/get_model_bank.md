# Get the session model bank

Retrieves the model bank set by
[`set_model_bank()`](https://sims1253.github.io/bayesim/reference/set_model_bank.md).
Returns NULL when no bank is active (e.g. `precompile = FALSE` or before
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
builds it).

## Usage

``` r
get_model_bank()
```

## Value

A named list of `brmsfit` prefit objects, or NULL.
