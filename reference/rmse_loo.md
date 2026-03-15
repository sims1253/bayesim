# Calculate PSIS-loo rmse for a given brms fit

Calculate PSIS-loo rmse for a given brms fit

## Usage

``` r
rmse_loo(fit, psis_object = NULL, return_object = FALSE, yrep = NULL, ...)
```

## Arguments

- fit:

  brms fit to calculate RMSE for

- psis_object:

  PSIS object for weights (optional)

- return_object:

  If TRUE, return a custom_loo_object

- yrep:

  Posterior predictive samples (optional)

- ...:

  Additional arguments to be passed to update() in case of reloo

## Value

`custom_loo_object` object with rmse acting as elpd, or a list with
rmse_loo and se_rmse_loo.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package
fit <- brms::brm(y ~ x, data = mydata)
rmse_loo(fit)
} # }
```
