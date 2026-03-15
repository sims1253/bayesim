# Full inverse forward sampling supported SBC

Full inverse forward sampling supported SBC

## Usage

``` r
ifs_SBC(
  fit,
  n_sims,
  ppred_data_gen,
  precon_sample = NULL,
  lb = NULL,
  ub = NULL,
  truncate = FALSE,
  post_warmup_samples = 1000,
  ...
)
```

## Arguments

- fit:

  brmsfit object that is used to generate SBC datasets and fit SBC
  models

- n_sims:

  Number of SBC simulations/datasets.

- ppred_data_gen:

  Function that generates the data for the newdata argument of
  posterior_predict to generate SBC datasets.

- precon_sample:

  The dataset that was used to precondition the fit.

- lb:

  Lower bound for the outcome. Is used to resample.

- ub:

  Upper bound for the outcome. Is used to resample.

- truncate:

  If TRUE, will truncate outcome values, if FALSE, will resample
  instead.

- post_warmup_samples:

  Is passed to the sbc model updates and used to make sure resampling
  only uses valid samples. Increase in case of numerical instability.

- ...:

  further parameters. Currently only passed to ppred_data_gen.

## Value

A tibble with SBC summaries and diagnostics.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires SBC and brms packages
} # }
```
