# Experimental grammar: primed-prior calibration

This declaration translates the weighted priming procedure in [case 2 of
*Primed Priors for Simulation-Based Validation of Bayesian
Models*](https://github.com/sims1253/implicit-priors/blob/cc59f18cc5162cf1b0005f7df4b89ff4f779bf15/case_2/sim_nonzero-betas/R/generator.R#L9-L45).
Building the plan executes no model fitting. The run at the end is
disabled. The grammar is an interface experiment, separate from the
supported
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
engine.

The procedure has two fitting stages. First, fit a model to weighted
priming data. Draw parameters from that fit and generate a new dataset.
Then fit the model to the new dataset **together with the same weighted
priming data**. That last step preserves the conditioning used to draw
the parameters. Dropping the priming rows changes the statistical
procedure.

``` r

source(system.file("experimental", "study-grammar.R", package = "bayesim"))
```

## Declare preparation and generation

We use one predictor and proper normal priors to keep the translation
readable. The source uses 250 predictors, an R2D2 prior and custom Stan
quantities. It also has separate zero-coefficient and
nonzero-coefficient variants. This example preserves weighted priming
and reuse of a generating fit; it does not reproduce those models or
their published results. [Source model and
design](https://github.com/sims1253/implicit-priors/blob/cc59f18cc5162cf1b0005f7df4b89ff4f779bf15/case_2/sim_nonzero-betas/_targets.R#L26-L44).

Each condition identifies a priming weight and an outer repetition.
Preparation runs once per condition. Inner replicates share its priming
data and generating fit. Making the outer repetition explicit avoids
treating all generated datasets as independent preparations.

``` r

prepare_priming <- function(condition, context) {
  x <- stats::rnorm(condition$m)
  priming <- data.frame(
    y = 0.5 + 0.1 * x + stats::rnorm(condition$m),
    x = x,
    w = condition$weight
  )
  prior <- c(
    brms::set_prior("normal(0, 2)", class = "b"),
    brms::set_prior("normal(0, 2)", class = "Intercept"),
    brms::set_prior("exponential(1)", class = "sigma")
  )
  fit <- brms::brm(
    brms::bf(y | weights(w) ~ x),
    data = priming, family = stats::gaussian(), prior = prior,
    backend = "cmdstanr", chains = 4, cores = 1,
    warmup = 1000, iter = 2000, seed = context$seed,
    refresh = 0
  )
  list(priming = priming, fit = fit)
}

generate_primed <- function(condition, context) {
  prepared <- context$prepared
  generating_draws <- posterior::as_draws_matrix(prepared$fit)
  draw_id <- context$replicate
  if (draw_id > nrow(generating_draws)) {
    stop("More inner replicates requested than available generating draws.")
  }
  newdata <- data.frame(
    y = 0, x = stats::rnorm(condition$n), w = 1
  )
  prediction <- brms::posterior_predict(
    prepared$fit, newdata = newdata, draw_ids = draw_id
  )
  newdata$y <- as.numeric(prediction[1, ])
  targets <- c("b_Intercept", "b_x", "sigma")
  list(
    train = rbind(newdata, prepared$priming),
    priming = prepared$priming,
    truth = stats::setNames(
      as.numeric(generating_draws[draw_id, targets]), targets
    ),
    generating_draw_id = draw_id,
    preparation_id = condition$condition_id
  )
}
```

The same draw index selects both the parameter truth and its predictive
data. No cycling occurs if the requested replicate count exceeds the
available draws. These are MCMC draws: this example does not establish
their independence or validate the finite-sample distribution of the
resulting calibration ranks.

## Declare refitting and artifacts

Refitting updates the prepared model using the combined dataset. It
therefore uses the original prior and both likelihood contributions.
Reusing the model also allows brms to reuse its compiled Stan program.

``` r

refit_primed <- function(data, context) {
  stats::update(
    context$prepared$fit, newdata = data$train,
    chains = context$settings$chains, cores = 1,
    warmup = context$settings$warmup, iter = context$settings$iter,
    seed = context$seed, refresh = 0
  )
}

extract_primed_draws <- function(fit, data, context) {
  posterior::subset_draws(
    posterior::as_draws_array(fit),
    variable = names(data$truth)
  )
}

extract_primed_diagnostics <- function(fit, data, context) {
  draws <- extract_primed_draws(fit, data, context)
  posterior::summarise_draws(draws)
}

primed_method <- study_method(
  fit = refit_primed,
  extract = list(
    draws = extract_primed_draws,
    diagnostics = extract_primed_diagnostics
  ),
  settings = list(chains = 4L, warmup = 1000L, iter = 2000L),
  version = "gaussian-refit-1"
)
```

The draw artifact preserves chains. Diagnostics can therefore use chain
information even though the rank calculation below pools draws.
Preparation fits are available during execution through
`context$prepared`. A durable run stores preparation separately from the
method’s `fit` artifact.

## Declare ranks and their aggregation

The rank measurement needs posterior draws and the generated-data
artifact, which contains truth and priming provenance. It returns one
row per target. Ties are randomized using the measurement’s RNG stream.
`max_rank` records the number of pooled posterior draws for that fit.

``` r

measure_primed_ranks <- function(artifacts, context) {
  draws <- posterior::as_draws_matrix(artifacts$draws)
  truth <- artifacts$data$truth
  ranks <- vapply(names(truth), function(target) {
    values <- draws[, target]
    if (!all(is.finite(values)) || !is.finite(truth[[target]])) {
      stop("Rank calculation requires finite draws and truth.")
    }
    below <- sum(values < truth[[target]])
    tied <- sum(values == truth[[target]])
    below + sample.int(tied + 1L, size = 1L) - 1L
  }, numeric(1))
  data.frame(
    target = names(truth), value = unname(ranks),
    max_rank = nrow(draws)
  )
}

summarize_primed_ranks <- function(results, artifacts, context) {
  if (!nrow(results)) return(data.frame())
  groups <- split(results, results$target)
  do.call(rbind, lapply(groups, function(rows) {
    data.frame(
      target = rows$target[1],
      value = mean(rows$value / rows$max_rank),
      n_ranked = nrow(rows),
      n_attempted = nrow(context$attempts)
    )
  }))
}
```

The grouped calculation reports mean normalized ranks and their
denominators. It is a descriptive check. It does not replace an SBC ECDF
test or account for autocorrelation in generating and refitting draws.
The source’s custom Gamma discrepancy is an across-replicate calculation
that could use the same grouped seam after its numerical behavior is
validated. [Source discrepancy
calculation](https://github.com/sims1253/implicit-priors/blob/cc59f18cc5162cf1b0005f7df4b89ff4f779bf15/case_3/R/utils.R#L347-L355).

## Build the plan

``` r

conditions <- expand.grid(
  weight = c(0, 0.01, 0.1), outer = 1:2,
  KEEP.OUT.ATTRS = FALSE
)
conditions$condition_id <- paste0(
  "weight-", conditions$weight, "-outer-", conditions$outer
)
conditions$m <- 100L
conditions$n <- 100L

primed_study <- study(
  "weighted-priming-gaussian",
  conditions = conditions,
  prepare = prepare_priming,
  generate = generate_primed,
  version = "priming-procedure-1"
)
primed_study <- with_methods(primed_study, brms = primed_method)
primed_study <- with_measure(
  primed_study, "rank", measure_primed_ranks,
  needs = c("draws", "data")
)
primed_study <- with_comparison(
  primed_study, "mean_normalized_rank", summarize_primed_ranks,
  by = c("condition_id", "method_id")
)
primed_study <- with_retention(
  primed_study, artifacts = c("data", "draws", "diagnostics")
)
plan_study(primed_study, replicates = 5L, seed = 671126974L)
#> weighted-priming-gaussian
#> 6 conditions x 5 replicates = 30 datasets
#> 1 methods; 30 fits
#> Prepare once per condition.
#>                  name       needs              stage
#>                  rank draws, data            per fit
#>  mean_normalized_rank             saved measurements
#> Retain: data, draws, diagnostics
#> Discard after use: 
#> Plan only: no data, models or cache records have been evaluated.
```

This plan contains six preparation conditions and thirty generated
datasets. Each dataset receives one refit. Planning validates the
declarations without calling brms, preparing data or inspecting a saved
run. Five inner replicates are enough to inspect the workflow shape, not
to establish calibration.

## Run and revisit measurements

The following chunk requires brms and a working CmdStan installation. It
runs six preparation fits and thirty refits, so it is excluded from
vignette builds.

``` r

primed_run <- run_study(
  primed_study, replicates = 5L, seed = 671126974L,
  path = "primed-priors-run", workers = 1L
)
primed_assessment <- assess_study(primed_run)
```

Measurement rows and attempted-fit records are always saved. The
declared retention also saves truth, the exact priming rows and their
weights, selected posterior draws, and diagnostics. A later interval or
rank measurement can use those artifacts in principle; it cannot recover
an omitted posterior variable. A new predictive measurement usually
needs a retained native fit or predictive draws and their data. Saving a
small rank table does not supply either.

`assess_study()` reruns comparisons of existing measurement rows. It
does not add a new per-fit measure to a saved run. Changing the
declaration and running it again must not be described as artifact reuse
until the runner provides an explicit reuse operation. A declared new
measure must state its artifact needs, and a reuse operation must reject
missing retained artifacts. Retaining `fit` would retain the refit. The
durable run stores each condition’s preparation separately so its
priming data and generating fit can be reused during execution.

The source also contains bodyfat resampling, Gamma response bounds,
latent variable masking and hand-written Stan quantities. Those remain
scientific code owned by their studies. This vignette makes no claim
that a generic generator can reproduce them or that the historical
`ifs_SBC()` calls are compatible with the current supported bayesim
interface. [Historical
call](https://github.com/sims1253/implicit-priors/blob/cc59f18cc5162cf1b0005f7df4b89ff4f779bf15/case_1/gamma_case_study.qmd#L172-L193),
[latent-data
transformation](https://github.com/sims1253/implicit-priors/blob/cc59f18cc5162cf1b0005f7df4b89ff4f779bf15/case_3/R/utils.R#L112-L128).
