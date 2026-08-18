# Brms Fitter

Fitter implementation for brms models. Extends the abstract Fitter class
with brms-specific configuration properties.

## Usage

``` r
BrmsFitter(
  name = "brms",
  supports_predictions = TRUE,
  supports_log_lik = TRUE,
  supports_loo = TRUE,
  supports_epred = TRUE,
  backend = "cmdstanr",
  chains = 4L,
  iter = 2000L,
  warmup = 1000L,
  thin = 1L,
  refresh = 0L,
  silent = 2L,
  cores = 1L,
  precompile = TRUE,
  allow_default_priors = FALSE,
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

- supports_epred:

  Logical indicating if posterior expectation predictions are supported
  (inherited)

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
  task. When precompiling, specify priors explicitly: some brms defaults
  are derived from the template dataset and would otherwise remain
  embedded in the reused compiled model.

- allow_default_priors:

  Logical, default FALSE. When FALSE, precompiled model banks reject
  model specs without an explicit prior with a fatal
  [`bayesim_config_error()`](https://sims1253.github.io/bayesim/reference/bayesim_config_error.md):
  brms derives data-dependent default priors from the template dataset,
  and they stay embedded in the compiled model that the bank reuses for
  every task (the whole study would silently be fit with the template's
  priors). Set TRUE to permit brms data-derived default priors to be
  embedded from the template data (rarely what you want; a notice is
  emitted once per run). Ignored when `precompile` is FALSE.

- stan_args:

  Named list of Stan/brms arguments passed through to the fit, e.g.
  `list(adapt_delta = 0.95, max_treedepth = 12, init = 0.1, threads = 2)`.
  NULL (default) uses brms/Stan defaults.

## Value

An S7 BrmsFitter object
