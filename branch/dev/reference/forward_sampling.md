# Simulate a new dataset using forward sampling.

Simulate a new dataset using forward sampling.

## Usage

``` r
forward_sampling(fit, i, n, ...)

# S3 method for class 'list'
forward_sampling(fit, i, n, ...)
```

## Arguments

- fit:

  An object of class brmsfit or a list of brmsfit objects.

- i:

  The index of a single posterior draw to simulate a dataset for. The
  index is passed to
  [`rstantools::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)'s
  "draw_ids" argument.

- n:

  Number of observations to simulate.

- ...:

  Potential additional arguments passed to
  [`brms::posterior_predict()`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.html).

## Value

A data.frame containing n observations for each variable in the fit.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- brms::brm(y ~ x, data = data.frame(y = rnorm(10), x = rnorm(10)))
forward_sampling(fit, i = 1, n = 100)
} # }
```
