# Validate a Fitter Object

Validates that a fitter object correctly implements the bayesim Fitter
interface. This function checks that the object is a valid S7 Fitter
class instance with all required properties and methods implemented.

## Usage

``` r
validate_fitter(
  fitter,
  smoke_test = FALSE,
  verbose = FALSE,
  data_bundle = NULL,
  fit_spec = NULL,
  task_ctx = NULL
)
```

## Arguments

- fitter:

  An S7 Fitter object to validate

- smoke_test:

  Logical, if TRUE run a quick fit test with sample data to verify that
  methods work correctly end-to-end

- verbose:

  Logical, if TRUE print progress messages during validation

- data_bundle:

  Optional representative data bundle for the conformance run. Supplying
  this is recommended for custom fitters whose data contract is not a
  `y ~ x` regression.

- fit_spec:

  Optional representative fit specification passed to
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
  during conformance testing.

- task_ctx:

  Optional task context passed to
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
  during conformance testing.

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

- `supports_epred` property exists and is logical

**Method Checks:**

- [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
  and
  [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)
  are implemented (the core contract)

- optional methods are required only when their `supports_*` capability
  is `TRUE`; unsupported methods have safe defaults

- [`fit_diagnostics()`](https://sims1253.github.io/bayesim/reference/fit_diagnostics.md)
  may use the default empty-list implementation

**Smoke Test (when smoke_test = TRUE):**

- Creates simple lm-like test data

- Calls
  [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
  and verifies `bayesim_fit_result` structure

- Calls
  [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)
  and verifies matrix with colnames

- If `supports_predictions`, calls
  [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
  and verifies output

- If `supports_log_lik`, calls
  [`log_lik_matrix()`](https://sims1253.github.io/bayesim/reference/log_lik_matrix.md)
  and verifies matrix output

- If `supports_epred`, calls
  [`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md)
  and verifies matrix output

- Calls
  [`fit_diagnostics()`](https://sims1253.github.io/bayesim/reference/fit_diagnostics.md)
  and verifies list output

## See also

[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[MockFitter](https://sims1253.github.io/bayesim/reference/MockFitter.md),
[`validate_metric()`](https://sims1253.github.io/bayesim/reference/validate_metric.md)

## Examples

``` r
# Validate a built-in fitter (basic check only)
validate_fitter(LinearRegressionFitter())

# Full validation with an end-to-end smoke test
validate_fitter(LinearRegressionFitter(), smoke_test = TRUE)

if (FALSE) { # \dontrun{
# Use in your own fitter's tests
my_fitter <- MyCustomFitter()
validate_fitter(my_fitter, smoke_test = TRUE, verbose = TRUE)
} # }
```
