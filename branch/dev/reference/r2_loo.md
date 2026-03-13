# Calculate PSIS-loo R² for a given brms fit

Calculate PSIS-loo R² for a given brms fit

## Usage

``` r
r2_loo(fit, psis_object = NULL, yrep = NULL, return_object = FALSE, ...)
```

## Arguments

- fit:

  brms fit to calculate R² for

- psis_object:

  PSIS object for weights (optional)

- yrep:

  Posterior predictive samples (optional)

- return_object:

  If TRUE, return a custom_loo_object

- ...:

  Additional arguments

## Value

`custom_loo_object` object with R² acting as elpd, or a list with r2_loo
and se_r2_loo.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package
fit <- brms::brm(y ~ x, data = mydata)
r2_loo(fit)
} # }
```
