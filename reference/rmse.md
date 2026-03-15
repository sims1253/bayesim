# Calculate the root-mean-squared-error for given y and yrep.

If psis-weights are supplied, the corresponding psis-rmse is returned.

## Usage

``` r
rmse(y, yrep, weights = NULL)
```

## Arguments

- y:

  Vector of observed values

- yrep:

  Matrix of predicted values (draws x observations)

- weights:

  PSIS weights (optional)

## Value

rmse for the given y and yrep vectors

## Examples

``` r
if (FALSE) { # \dontrun{
y <- c(1, 2, 3)
yrep <- matrix(rnorm(300), nrow = 100, ncol = 3)
rmse(y, yrep)
} # }
```
