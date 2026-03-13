# Validate Fitter Interface

Validates that a fitter object conforms to the Fitter S7 class
interface. The fitter must be an S7 object that inherits from the Fitter
class and implements all required methods.

## Usage

``` r
validate_fitter_interface(fitter)
```

## Arguments

- fitter:

  An S7 object to validate as a Fitter.

## Value

The input `fitter`, invisibly, if validation passes.

## Details

The fitter must satisfy the following requirements:

- Must be an S7 object (checked via S7::S7_inherits())

- Must inherit from the "Fitter" class

- Must implement the following methods:

  - `fit(data_bundle, fit_spec, seed, task_ctx)`

  - `extract_draws(fit_result, variables = NULL)`

  - `predict(fit_result, newdata = NULL, seed = NULL)`

  - `log_lik(fit_result, newdata = NULL)`

  - `loo(fit_result)`

  - `diagnostics(fit_result)`

## Errors

Throws a `bayesim_contract_error` condition if validation fails.

## See also

[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[MockFitter](https://sims1253.github.io/bayesim/reference/MockFitter.md)

## Examples

``` r
# Validate the mock fitter
mock_fitter <- MockFitter()
validate_fitter_interface(mock_fitter)
```
