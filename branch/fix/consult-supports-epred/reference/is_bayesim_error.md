# Check if a condition is a bayesim error

Check if a condition is a bayesim error

## Usage

``` r
is_bayesim_error(cond)
```

## Arguments

- cond:

  A condition object to test

## Value

TRUE if the condition is a bayesim error, FALSE otherwise

## Examples

``` r
if (FALSE) { # \dontrun{
tryCatch(
  bayesim_config_error("Invalid config"),
  bayesim_error = function(cond) {
    is_bayesim_error(cond) # TRUE
  }
)
} # }
```
