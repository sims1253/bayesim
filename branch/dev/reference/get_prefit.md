# Create a precompiled brms model object

Creates a brmsfit object with compiled Stan code that can be efficiently
updated with new data during simulation runs, avoiding recompilation.

## Usage

``` r
get_prefit(fit_conf, stan_pars)
```

## Arguments

- fit_conf:

  A list or data.frame row containing model configuration:

  - fit_family: The model family (e.g., "gaussian", "student_t")

  - fit_link: The link function (e.g., "identity", "log")

  - formula: The model formula as a character string

  - prior: (optional) Prior specifications

- stan_pars:

  A list containing Stan parameters, must include:

  - backend: Stan backend ("cmdstanr" or "rstan")

## Value

A brmsfit object with compiled Stan code (chains = 0)

## Examples

``` r
if (FALSE) { # \dontrun{
fit_conf <- list(
  fit_family = "gaussian",
  fit_link = "identity",
  formula = "y ~ x"
)
stan_pars <- list(backend = "cmdstanr")
prefit <- get_prefit(fit_conf, stan_pars)
} # }
```
