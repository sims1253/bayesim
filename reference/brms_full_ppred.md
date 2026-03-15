# Full forward sampling of a the response of brms fit, including multivariate models.

Full forward sampling of a the response of brms fit, including
multivariate models.

## Usage

``` r
brms_full_ppred(fit, newdata = NULL, draws = NULL, validate_all = FALSE)
```

## Source

This function is taken from the SBC package
(https://github.com/hyunjimoon/SBC) and only here to ensure stability,
as it is not exported in SBC.

## Arguments

- fit:

  An object of class `brmsfit`

- newdata:

  An optional data.frame for which to evaluate predictions. If NULL
  (default), the original data of the model is used.

- draws:

  An integer vector specifying the posterior draws to be used. If NULL
  (the default), all draws are used.

- validate_all:

  if TRUE, validation of input data will be done in all iterations,
  otherwise only once

## Value

A list of data.frames containing the draws.

## Examples

``` r
# Pending
```
