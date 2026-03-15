# Calculate R² on test data

Calculate R² on test data

## Usage

``` r
r2_test(fit, newdata, return_object = FALSE)
```

## Arguments

- fit:

  A brmsfit object

- newdata:

  Test data frame

- return_object:

  If TRUE, return a custom_loo_object

## Value

A named list with r2_test and se_r2_test, or a custom_loo_object

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package
fit <- brms::brm(y ~ x, data = train_data)
r2_test(fit, test_data)
} # }
```
