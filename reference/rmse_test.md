# Calculate RMSE on test data

Calculate RMSE on test data

## Usage

``` r
rmse_test(fit, newdata, return_object = FALSE)
```

## Arguments

- fit:

  A brmsfit object

- newdata:

  Test data frame

- return_object:

  If TRUE, return a custom_loo_object

## Value

A named list with rmse_test and se_rmse_test, or a custom_loo_object

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package
fit <- brms::brm(y ~ x, data = train_data)
rmse_test(fit, test_data)
} # }
```
