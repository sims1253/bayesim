# Validate Fit Result Interface

Validates that a fit_result conforms to the bayesim_fit_result
interface. This is a wrapper around
[`validate_bayesim_fit_result()`](https://sims1253.github.io/bayesim/reference/validate_bayesim_fit_result.md)
that provides clear error messages for contract violations.

## Usage

``` r
validate_fit_result_interface(fit_result)
```

## Arguments

- fit_result:

  A bayesim_fit_result object to validate.

## Value

The input `fit_result`, invisibly, if validation passes.

## Details

The fit_result must satisfy the bayesim_fit_result class requirements:

- Must inherit from "bayesim_fit_result" class

- If `success` is TRUE, `error` must be NULL

- If `success` is FALSE, `error` must be non-NULL

- `timing$total` must be a non-negative scalar numeric

- If `draws` is not NULL, it must be a matrix with column names

- `warnings` must be a character vector

- `diagnostics` must be a list

## Errors

Throws a `bayesim_contract_error` condition if validation fails.

## See also

[`validate_bayesim_fit_result()`](https://sims1253.github.io/bayesim/reference/validate_bayesim_fit_result.md),
[`new_fit_result()`](https://sims1253.github.io/bayesim/reference/new_fit_result.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a valid fit result
draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
colnames(draws) <- c("alpha", "beta")
result <- new_fit_result(
  success = TRUE,
  draws = draws,
  diagnostics = list(rhat_max = 1.01),
  timing = list(total = 5.0, warmup = 2.5, sample = 2.5)
)
validate_fit_result_interface(result)
} # }
```
