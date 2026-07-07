# CmdStan Fitter (user-supplied Stan programs)

Run user-supplied Stan programs via cmdstanr without brms. Compilation
is cached by cmdstanr (by file hash); each daemon compiles-or-cache-hits
on first use, so there is no model-bank integration (D2).

The Stan program may declare generated-quantities blocks for `log_lik`
(a vector of pointwise log-likelihoods) and optionally `epred` (the
expectation `mu`). Declare their names via the `log_lik` / `epred`
arguments.

**Newdata prediction is out of scope** (raw Stan has no newdata
semantics): `supports_predictions` is FALSE unless `epred` is given, in
which case `predict_epred` returns the in-sample GQ matrix and
`predict_fit` is unsupported. Test-set metrics require brms or a custom
fitter.

## Usage

``` r
CmdStanFitter(
  stan_file = NULL,
  stan_code = NULL,
  stan_data,
  log_lik = NULL,
  epred = NULL,
  chains = 4L,
  iter_warmup = 1000L,
  iter_sampling = 1000L,
  adapt_delta = NULL,
  max_treedepth = NULL,
  parallel_chains = 1L,
  init = NULL,
  ...
)
```

## Arguments

- stan_file:

  Path to a `.stan` file, or NULL if `stan_code` is supplied.

- stan_code:

  Character string of Stan code (used when `stan_file` is NULL). Either
  `stan_file` or `stan_code` is required.

- stan_data:

  Function `function(data_bundle, fit_spec) -> list` mapping a bayesim
  data bundle + fit spec to the Stan `data` list. Required.

- log_lik:

  Name of the generated-quantities log-lik vector in the Stan program,
  or NULL if the program has no log_lik GQ. When NULL,
  `supports_log_lik`/`supports_loo` are FALSE and elpd metrics degrade.

- epred:

  Optional name of an epred/mu GQ matrix or vector. When supplied,
  `supports_predictions` is TRUE and `predict_epred` returns the
  in-sample GQ matrix (S x N).

- chains:

  Integer number of MCMC chains.

- iter_warmup:

  Integer warmup iterations per chain.

- iter_sampling:

  Integer sampling iterations per chain.

- adapt_delta:

  Numeric adapt_delta control parameter, or NULL for Stan default.

- max_treedepth:

  Integer max_treedepth control parameter, or NULL.

- parallel_chains:

  Integer number of chains to run in parallel.

- init:

  Passed through to `$sample()` (NULL, "random", "0", or a list).

- ...:

  Additional named arguments passed through to cmdstanr's `$sample()`.

## Value

An S7 `CmdStanFitter` object.

## See also

[Fitter](https://sims1253.github.io/bayesim/reference/Fitter.md),
[BrmsFitter](https://sims1253.github.io/bayesim/reference/BrmsFitter.md),
[LinearRegressionFitter](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fitter <- CmdStanFitter(
  stan_file = "model.stan",
  stan_data = function(bundle, spec) list(N = nrow(bundle$train), y = bundle$train$y),
  log_lik = "log_lik"
)
} # }
```
