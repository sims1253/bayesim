# Prepare a list of brmsfit objects to update during repeated simulations

This is mainly used to save on constant recompilation times.

## Usage

``` r
build_prefit_list(fit_configuration, stan_pars)
```

## Arguments

- fit_configuration:

  A named list that currently holds family, link and prior.

- stan_pars:

  A named list which contains a backend field.

## Value

A named list of precompiled fit objects keyed by fit configuration.
