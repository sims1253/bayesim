# Generate lookup keys for fit configurations to retrieve prefit objects matching said config.

Generate lookup keys for fit configurations to retrieve prefit objects
matching said config.

## Usage

``` r
fit_conf_key(fit_conf)
```

## Arguments

- fit_conf:

  A list containing fit configuration with elements: fit_family,
  fit_link, formula, and prior

## Value

A hash generated from the fit configuration

## Examples

``` r
fit_conf_key(
  list(
    fit_family = "gaussian",
    fit_link = "identity",
    prior = list(c(brms::set_prior("", class = "Intercept")))
  )
)
#> [1] "45a2f18e5dd478acd6feea29f6d3a3b6"
```
