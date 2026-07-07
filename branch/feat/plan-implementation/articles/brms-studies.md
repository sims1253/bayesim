# brms studies: the model bank and model_grid()

``` r

library(bayesim)
```

This vignette covers running brms-backed simulation studies: how the
model bank eliminates per-task compilation, how to declare
model-comparison grids with
[`model_grid()`](https://sims1253.github.io/bayesim/reference/model_grid.md)
/
[`brms_model()`](https://sims1253.github.io/bayesim/reference/brms_model.md),
how to pass sampler arguments, and how to keep memory flat on long runs.
All chunks require brms and a working CmdStan install, so they are shown
but not evaluated here.

## The model bank: compile once, fit thousands of times

A naive brms simulation study recompiles the Stan model for every task —
minutes of C++ compilation per seconds of sampling.
[`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
avoids this with a *model bank*: at the start of
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md),
each distinct model spec in the `fit_grid` is compiled once (via
`brms::brm(chains = 0)`), and every task reuses the compiled binary
through `stats::update(recompile = FALSE)`.

Two properties control this:

- `precompile` (default `TRUE`): build the bank at run start. Set to
  `FALSE` to fall back to a fresh
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) per
  task (only sensible for tiny studies or debugging).
- The bank is shipped to mirai daemons once per run, so parallel workers
  reuse the same binaries.

brms does **not** warn when `recompile = FALSE` is used against
structurally incompatible data — it would silently reuse the binary
against the wrong model frame. bayesim therefore compares the Stan data
*structure* (via
[`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html))
between the compiled template and each task’s data, and aborts the run
loudly on a mismatch. If your `data_grid` rows produce data with
different shapes (e.g. varying factor levels), either make them
structurally identical or set `precompile = FALSE`.

## Declaring models: `brms_model()` and `model_grid()`

Hand-building list-columns (`fit_grid$formula <- list(...)`) works but
is error-prone.
[`model_grid()`](https://sims1253.github.io/bayesim/reference/model_grid.md)
assembles a tidy fit grid from named
[`brms_model()`](https://sims1253.github.io/bayesim/reference/brms_model.md)
specs, validating each at construction so mistakes surface before any
compilation:

``` r

fit_grid <- model_grid(
  gaussian  = brms_model(y ~ x, brms::brmsfamily("gaussian")),
  student   = brms_model(y ~ x, brms::brmsfamily("student")),
  lognormal = brms_model(y ~ x, brms::brmsfamily("lognormal"))
)
fit_grid
#> # A tibble: 3 x 6
#>   model     formula  family  prior  stanvars stan_file
#>   <chr>     <list>   <list>  <list> <list>   <list>
#> 1 gaussian  <formula> <family> <NULL> <NULL>  <NULL>
#> 2 student   <formula> <family> <NULL> <NULL>  <NULL>
#> 3 lognormal <formula> <family> <NULL> <NULL>  <NULL>
```

The `model` name column lands in the result summary as `fit_model`, so
per-model aggregation is `summarize_simulation(result)` or
`performance_measures(result)` with no extra arguments.

## Case study: comparing likelihoods for skewed data

Which likelihood best describes right-skewed data — Gaussian, Student-t,
or lognormal? Generate skewed data, fit all three models, and compare
expected log predictive density (ELPD):

``` r

skew_generator <- function(data_spec, task_ctx) {
  n <- data_spec$n
  x <- stats::rnorm(n)
  mu <- data_spec$beta * x
  y <- mu + stats::rgamma(n, shape = 2, scale = 1)  # right-skewed noise
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(beta = data_spec$beta),
    vars_of_interest = "beta"
  )
}

config <- simulation_config(
  data_grid = data.frame(n = 200, beta = 1),
  fit_grid = model_grid(
    gaussian  = brms_model(y ~ x, brms::brmsfamily("gaussian")),
    student   = brms_model(y ~ x, brms::brmsfamily("student")),
    lognormal = brms_model(y ~ x, brms::brmsfamily("lognormal"))
  ),
  data_generator = skew_generator,
  fitter = BrmsFitter(chains = 2L, iter = 1000L, warmup = 500L),
  metrics = list(
    elpd_loo_metric(),
    posterior_summary_metric(),
    convergence_metric()
  ),
  n_replicates = 50L,
  seed = 42L
)

result <- run_simulation(config, workers = 4)
summarize_simulation(result, metrics = "elpd_loo__elpd")
```

Three models compile exactly once each; all 150 fits reuse the binaries.

## Sampler arguments: `stan_args`

`BrmsFitter(stan_args = ...)` passes sampler controls through to
brms/Stan. `adapt_delta` and `max_treedepth` are mapped into the
`control` list; `init` and `threads` pass through directly:

``` r

fitter <- BrmsFitter(
  chains = 4L,
  stan_args = list(adapt_delta = 0.95, max_treedepth = 12, init = 0.1)
)
```

## Warning-conditional retention

Diagnosing occasional divergent transitions across thousands of tasks
needs the fit object — but only for the tasks that misbehaved. Retention
specs can be conditional on warnings, keeping heavy artifacts only where
something went wrong:

``` r

config <- simulation_config(
  ...,
  retain = list(
    success = c("metrics", "diagnostics"),
    warning = c("metrics", "diagnostics", "fit", "draws")
  )
)
```

Clean tasks stay light; tasks that emitted sampler warnings keep their
full fit for post-hoc inspection.

## SBC with brms generators

For calibration checking of brms models,
[`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md)
draws the truth from the model prior (a `sample_prior = "only"` fit) and
[`ifs_generator()`](https://sims1253.github.io/bayesim/reference/ifs_generator.md)
draws it from a preconditioning posterior. Both forward-simulate the
response — including dependency-ordered simulation for multivariate
models. See
[`vignette("sbc-and-calibration")`](https://sims1253.github.io/bayesim/articles/sbc-and-calibration.md)
for the full workflow.

## Next steps

- [`vignette("parallel-and-hpc")`](https://sims1253.github.io/bayesim/articles/parallel-and-hpc.md)
  — daemons, checkpoint/resume, memory.
- [`vignette("custom-fitters")`](https://sims1253.github.io/bayesim/articles/custom-fitters.md)
  — raw Stan via
  [`CmdStanFitter()`](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md)
  when brms cannot express your model.
- [`vignette("reproducibility")`](https://sims1253.github.io/bayesim/articles/reproducibility.md)
  — what the config fingerprint covers.
