# Compute ELPD on test data

Compute ELPD on test data

## Usage

``` r
elpd_test(fit, newdata, return_object = FALSE)
```

## Arguments

- fit:

  A brmsfit object

- newdata:

  Test data frame

- return_object:

  If TRUE, return a custom_loo_object

## Value

A named list with elpd_test and se_elpd_test, or a custom_loo_object

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package
fit <- brms::brm(y ~ x, data = train_data)
elpd_test(fit, test_data)
} # }
```
