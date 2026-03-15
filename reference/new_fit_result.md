# Create a new bayesim_fit_result object

Constructs a result object from a Bayesian model fitting operation.

## Usage

``` r
new_fit_result(
  success = TRUE,
  fit = NULL,
  draws = NULL,
  diagnostics = list(),
  timing = list(total = 0, warmup = 0, sample = 0),
  warnings = character(),
  error = NULL,
  data_bundle = NULL
)
```

## Arguments

- success:

  Logical scalar indicating if the fit succeeded

- fit:

  The backend fit object (may be NULL if removed by retention policy)

- draws:

  Matrix of posterior draws (S x P), with column names for parameters

- diagnostics:

  Named list of diagnostic values (scalars or named numeric vectors)

- timing:

  List containing total, warmup, and sample timing in seconds

- warnings:

  Character vector of warning messages captured during fitting

- error:

  A condition object if the fit failed, NULL otherwise

- data_bundle:

  Optional list containing the training/test data used for fitting

## Value

A validated `bayesim_fit_result` object

## Details

The bayesim_fit_result class encapsulates all outputs from a Bayesian
model fitting operation, including the posterior draws, diagnostics,
timing information, and any warnings or errors encountered.

Validation rules:

- If `success` is FALSE, `error` must be non-NULL

- If `success` is TRUE, `error` must be NULL

- `timing$total` must be non-negative

- If `draws` is not NULL, it must be a matrix with column names

## Examples

``` r
# Successful fit
draws <- matrix(rnorm(1000), ncol = 2, nrow = 500)
colnames(draws) <- c("alpha", "beta")
result <- new_fit_result(
  success = TRUE,
  draws = draws,
  diagnostics = list(rhat = c(alpha = 1.01, beta = 1.00)),
  timing = list(total = 10.5, warmup = 5.0, sample = 5.5)
)

# Failed fit
result <- new_fit_result(
  success = FALSE,
  error = simpleError("Convergence failed"),
  diagnostics = list(),
  timing = list(total = 2.0, warmup = 2.0, sample = 0)
)
```
