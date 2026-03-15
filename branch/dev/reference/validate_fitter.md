# Validate a Fitter Object

Validates that a fitter object correctly implements the bayesim Fitter
interface. This function checks that the object is a valid S7 Fitter
class instance with all required properties and methods implemented.

## Usage

``` r
validate_fitter(fitter, smoke_test = FALSE, verbose = FALSE)
```

## Arguments

- fitter:

  An S7 Fitter object to validate

- smoke_test:

  Logical, if TRUE run a quick fit test with sample data to verify that
  methods work correctly end-to-end

- verbose:

  Logical, if TRUE print progress messages during validation

## Value

The validated fitter object (invisibly) if valid, otherwise raises an
error with details about what failed

## Details

The validation performs the following checks:

**Property Checks:**

- Object is an S7 Fitter class

- `name` property exists and is character

- `supports_predictions` property exists and is logical

- `supports_log_lik` property exists and is logical

- `supports_loo` property exists and is logical

**Method Checks:**

- [`fit()`](https://sims1253.github.io/bayesim/reference/fit.md) method
  is implemented

- [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)
  method is implemented

- [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
  method is implemented

- [`log_lik()`](https://sims1253.github.io/bayesim/reference/log_lik.md)
  method is implemented

- [`loo()`](https://sims1253.github.io/bayesim/reference/loo.md) method
  is implemented

- [`diagnostics()`](https://sims1253.github.io/bayesim/reference/diagnostics.md)
  method is implemented

**Smoke Test (when smoke_test = TRUE):**

- Creates simple lm-like test data

- Calls [`fit()`](https://sims1253.github.io/bayesim/reference/fit.md)
  and verifies `bayesim_fit_result` structure

- Calls
  [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)
  and verifies matrix with colnames

- If `supports_predictions`, calls
  [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
  and verifies output

- If `supports_log_lik`, calls
  [`log_lik()`](https://sims1253.github.io/bayesim/reference/log_lik.md)
  and verifies matrix output

- Calls
  [`diagnostics()`](https://sims1253.github.io/bayesim/reference/diagnostics.md)
  and verifies list output

## See also

[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[MockFitter](https://sims1253.github.io/bayesim/reference/MockFitter.md)

## Examples

``` r
# Validate the built-in MockFitter (basic check only)
validate_fitter(MockFitter())

# Full validation with smoke test
validate_fitter(MockFitter(), smoke_test = TRUE, verbose = TRUE)
#> Checking S7 Fitter class...
#>   [OK] Object is an S7 Fitter class
#> Checking required properties...
#>   [OK] Property 'name' exists: mock
#>   [OK] Property 'supports_predictions' exists: TRUE
#>   [OK] Property 'supports_log_lik' exists: TRUE
#>   [OK] Property 'supports_loo' exists: TRUE
#> Checking required S7 methods...
#>   [OK] Method 'fit()' is implemented
#>   [OK] Method 'extract_draws()' is implemented
#>   [OK] Method 'predict_fit()' is implemented
#>   [OK] Method 'log_lik()' is implemented
#>   [OK] Method 'loo()' is implemented
#>   [OK] Method 'diagnostics()' is implemented
#> Running smoke test with sample data...
#>   Testing fit()...
#>     [OK] fit() returns valid bayesim_fit_result
#>   Testing extract_draws()...
#>     [OK] extract_draws() returns matrix with colnames
#>   Testing predict_fit()...
#>     [OK] predict_fit() returns valid predictions
#>   Testing log_lik()...
#>     [OK] log_lik() returns valid matrix
#>   Testing diagnostics()...
#>     [OK] diagnostics() returns list
#> Smoke test completed successfully!
#> Validation passed!

# Use in your fitter tests
if (FALSE) { # \dontrun{
my_fitter <- MyCustomFitter()
validate_fitter(my_fitter, smoke_test = TRUE)
cat("Fitter is valid!\n")
} # }
```
