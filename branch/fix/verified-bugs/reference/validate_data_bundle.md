# Validate a Data Bundle

Validates that a data_bundle conforms to the required structure for use
in simulation tasks. The data_bundle is the output of a data_generator
function and contains all data-related objects needed for model fitting
and metric computation.

## Usage

``` r
validate_data_bundle(data_bundle)
```

## Arguments

- data_bundle:

  A list containing data and metadata for a simulation task.

## Value

The input `data_bundle`, invisibly, if validation passes.

## Details

The data_bundle must have the following structure:

- `train`: A data.frame with at least 1 row (required)

- `test`: NULL or a data.frame (optional)

- `response`: A scalar character naming the response column in train
  (and test if not NULL)

- `true_params`: A named numeric vector where names exactly match
  `vars_of_interest`

- `vars_of_interest`: A non-empty character vector of unique names

- `meta`: Optional named list with scalar values only

Validation rules:

- `train` must be a data.frame with nrow \>= 1

- `test` must be NULL or a data.frame

- `response` must be a scalar character present in train (and test if
  not NULL)

- `true_params` must be a named numeric vector

- `vars_of_interest` must be a non-empty unique character vector

- `setequal(names(true_params), vars_of_interest)` must be TRUE

- No duplicate names in any named vector/list

- `meta` must be a named list with scalar values only (if present)

## Errors

Throws a `bayesim_data_error` condition if validation fails.

## Examples

``` r
# Valid data bundle
data_bundle <- list(
  train = data.frame(x = 1:10, y = rnorm(10)),
  test = data.frame(x = 11:15, y = rnorm(5)),
  response = "y",
  true_params = c(beta = 1.5, sigma = 0.5),
  vars_of_interest = c("beta", "sigma")
)
validate_data_bundle(data_bundle)
```
