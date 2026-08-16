# Extract the response variable name from a brmsfit. Falls back through formula forms since sample_prior="only" fits may strip the standard formula() variables.

Extract the response variable name from a brmsfit. Falls back through
formula forms since sample_prior="only" fits may strip the standard
formula() variables.

## Usage

``` r
.fit_response_name(fit)
```
