# Extract Fit Diagnostics

Extract convergence and fit diagnostics.

## Usage

``` r
diagnostics(fitter, fit_result)
```

## Arguments

- fitter:

  An S7 Fitter object

- fit_result:

  A `bayesim_fit_result` object from
  [`fit()`](https://sims1253.github.io/bayesim/reference/fit.md)

## Value

A named list of scalar diagnostic values, which may include:

- `rhat_max`: Maximum R-hat statistic across parameters

- `ess_bulk`: Minimum bulk effective sample size

- `ess_tail`: Minimum tail effective sample size

- `divergent`: Number of divergent transitions

- `max_treedepth`: Number of max treedepth warnings

- Fitter-specific diagnostics
