# Brms Fitter

Fitter implementation for brms models.

## Usage

``` r
BrmsFitter(
  name = character(0),
  supports_predictions = TRUE,
  supports_log_lik = TRUE,
  supports_loo = TRUE,
  backend = "cmdstanr",
  chains = 4L,
  iter = 2000L,
  warmup = 1000L,
  thin = 1L,
  refresh = 0L,
  silent = 2L,
  cores = 1L
)
```

## Arguments

- name:

  Character string identifying the fitter (inherited from Fitter)

- supports_predictions:

  Logical indicating if predictions are supported (inherited)

- supports_log_lik:

  Logical indicating if log-likelihood is supported (inherited)

- supports_loo:

  Logical indicating if LOO-CV is supported (inherited)

- backend:

  Character string for Stan backend ("cmdstanr" or "rstan")

- chains:

  Integer number of MCMC chains

- iter:

  Integer total iterations per chain

- warmup:

  Integer warmup iterations per chain

- thin:

  Integer thinning interval

- refresh:

  Integer refresh rate for progress output

- silent:

  Integer verbosity level (0, 1, or 2)

- cores:

  Integer number of cores for parallel processing

## Value

An S7 BrmsFitter object
