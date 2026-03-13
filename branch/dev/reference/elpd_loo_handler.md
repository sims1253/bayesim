# Extract ELPD-LOO from a brms fit

Extract ELPD-LOO from a brms fit

## Usage

``` r
elpd_loo_handler(fit)
```

## Arguments

- fit:

  A brmsfit object

## Value

A named list with p_loo, se_p_loo, elpd_loo, se_elpd_loo, looic,
se_looic

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package and a fitted model
fit <- brms::brm(y ~ x, data = mydata)
elpd_loo_handler(fit)
} # }
```
