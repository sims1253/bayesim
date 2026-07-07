# Brms Fitter

Fitter implementation for brms models. Extends the abstract Fitter class
with brms-specific configuration properties.

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
  cores = 1L,
  precompile = TRUE,
  stan_args = list()
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

- precompile:

  Logical; if TRUE (default), the model bank compiles each distinct
  model spec once via `brms::brm(chains = 0)` and reuses the prefit
  across tasks via `stats::update(recompile = FALSE)`. Set to FALSE to
  fall back to a fresh
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) per
  task.

- stan_args:

  Named list of Stan/brms arguments passed through to the fit, e.g.
  `list(adapt_delta = 0.95, max_treedepth = 12, init = 0.1, threads = 2)`.
  NULL (default) uses brms/Stan defaults.

## Value

An S7 BrmsFitter object
