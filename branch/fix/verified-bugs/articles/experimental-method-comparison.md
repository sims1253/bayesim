# Experimental study grammar: compare prior specifications

This example declares a small regression study with two shrinkage priors
and one stricter sampler setting. It uses the direct Stan programs from
[GDR2 at commit
`3b1350a`](https://github.com/jear2412/GDR2/tree/3b1350a7d56dbc67af6e31d7ea7ce114ae72dba3).
The design is an illustration of the experimental grammar, not a
reproduction of the paper’s simulations. Building this vignette declares
and plans the study; it does not compile models or sample posteriors.

The GDR2 example compares prior specifications using a shared regression
[training and test
dataset](https://github.com/jear2412/GDR2/blob/3b1350a7d56dbc67af6e31d7ea7ce114ae72dba3/sim/example.Rmd#L57-L167).
Here the likelihood stays Gaussian. The method IDs distinguish both the
prior and sampler setting, so either can vary without changing the data
generator.

``` r

library(bayesim)
source(system.file("experimental", "study-grammar.R", package = "bayesim"))
```

## Data and truth

Each replicate contains a training set, a held-out test set, and named
coefficient truths. Conditions vary the training sample size. All
methods receive the same generated object for a given `dataset_id`.

The nonzero coefficients are fixed across replicates. This is a recovery
and prediction study; uniform SBC ranks are not an expected result of
this design.

``` r

generate_regression <- function(condition, context) {
  p <- condition$p
  beta <- c(1, -0.5, rep(0, p - 2L))
  names(beta) <- sprintf("beta[%d]", seq_len(p))
  intercept <- 0.5
  sigma <- 1
  make_sample <- function(n) {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    y <- as.numeric(intercept + x %*% beta + stats::rnorm(n, sd = sigma))
    list(x = x, y = y)
  }
  list(
    train = make_sample(condition$n),
    test = make_sample(condition$n_test),
    truth = c(beta, b_Intercept = intercept, sigma = sigma)
  )
}

conditions <- data.frame(
  condition_id = c("n40", "n100"),
  n = c(40L, 100L),
  p = 6L,
  n_test = 100L
)
```

## Direct Stan methods

Obtain `stan/models/R2D2.stan` and `stan/models/logitR2.stan` from the
pinned GDR2 checkout before running the study. Set `BAYESIM_GDR2_ROOT`
to its root; the default path below is `GDR2` in the working directory.
Planning requires neither these files nor CmdStan. Execution requires
both files, cmdstanr, posterior, and a working CmdStan installation. No
code below downloads models.

Both source programs include the intercept in `p` and `X`, but their
`beta` vectors contain slopes only. The data mapping follows the
[original
helper](https://github.com/jear2412/GDR2/blob/3b1350a7d56dbc67af6e31d7ea7ce114ae72dba3/utils/utils.R#L23-L51).
The programs center the training response and expose the intercept on
the original scale as `b_Intercept`.

``` r

fit_gdr2 <- function(data, context) {
  settings <- context$settings
  if (!file.exists(settings$stan_file)) {
    stop("Missing GDR2 Stan source: ", settings$stan_file)
  }
  p <- ncol(data$train$x)
  stan_data <- list(
    N = length(data$train$y),
    p = p + 1L,
    X = cbind(1, data$train$x),
    y = data$train$y,
    Ntest = length(data$test$y),
    Xtest = cbind(1, data$test$x),
    ytest = data$test$y,
    scale_sigma = stats::sd(data$train$y),
    prior_only = 0L,
    R2_mean = settings$R2_mean,
    R2_prec = settings$R2_prec
  )
  if (settings$prior == "dirichlet") {
    stan_data$R2_alpha <- rep(settings$concentration, p)
  } else if (settings$prior == "logit_normal") {
    stan_data$mu_logitphi <- rep(0, p)
    # Stan expects a lower Cholesky factor; a diagonal scale is valid.
    stan_data$sigma_logitphi <- diag(settings$logit_sd, p)
  } else {
    stop("Unknown prior specification: ", settings$prior)
  }
  model <- cmdstanr::cmdstan_model(settings$stan_file, force_recompile = FALSE)
  fit <- model$sample(
    data = stan_data,
    seed = context$seed,
    chains = settings$chains,
    parallel_chains = 1L,
    iter_warmup = settings$warmup,
    iter_sampling = settings$sampling,
    adapt_delta = settings$adapt_delta,
    refresh = 0
  )
  if (any(fit$return_codes() != 0L)) stop("A CmdStan chain failed")
  fit
}
```

Compilation and executable reuse follow cmdstanr’s file handling. The
study version below records which external source is intended; the
example does not verify a checkout or hash its files. Change the version
when changing that source. One chain runs at a time within each method
to keep outer study workers from each launching several chain processes.

Extract only the artifacts the measurements need. Posterior draws keep
their iteration and chain dimensions. Test scores carry observation IDs
so the comparison can check that the methods scored the same rows. The
stable log-mean-exp calculation averages predictive densities over
posterior draws.

``` r

extract_parameters <- function(fit, data, context) {
  fit$draws(variables = names(data$truth), format = "draws_array")
}

extract_test_scores <- function(fit, data, context) {
  variables <- sprintf("log_lik_test[%d]", seq_along(data$test$y))
  draws <- fit$draws(variables = variables, format = "draws_matrix")
  draws <- draws[, variables, drop = FALSE]
  log_mean_exp <- function(x) {
    if (anyNA(x) || any(x == Inf)) stop("Invalid predictive log density")
    if (all(x == -Inf)) return(-Inf)
    center <- max(x)
    center + log(mean(exp(x - center)))
  }
  data.frame(
    observation_id = seq_along(data$test$y),
    log_score = apply(draws, 2L, log_mean_exp)
  )
}

extract_diagnostics <- function(fit, data, context) {
  fit$diagnostic_summary()
}

make_gdr2_method <- function(stan_file, prior, adapt_delta = 0.90) {
  study_method(
    fit = fit_gdr2,
    extract = list(
      draws = extract_parameters,
      test_scores = extract_test_scores,
      diagnostics = extract_diagnostics
    ),
    settings = list(
      stan_file = stan_file, prior = prior,
      R2_mean = 0.5, R2_prec = 1,
      concentration = 0.5, logit_sd = 1,
      chains = 4L, warmup = 1000L, sampling = 1000L,
      adapt_delta = adapt_delta
    ),
    version = "GDR2-3b1350a7d56dbc67af6e31d7ea7ce114ae72dba3-v1"
  )
}

model_root <- Sys.getenv("BAYESIM_GDR2_ROOT", unset = "GDR2")
method_study <- study(
  "gdr2-prior-and-sampler",
  conditions = conditions,
  generate = generate_regression
)
method_study <- with_methods(
  method_study,
  dirichlet = make_gdr2_method(
    file.path(model_root, "stan/models/R2D2.stan"), "dirichlet"
  ),
  logit_normal = make_gdr2_method(
    file.path(model_root, "stan/models/logitR2.stan"), "logit_normal"
  ),
  logit_normal_strict = make_gdr2_method(
    file.path(model_root, "stan/models/logitR2.stan"), "logit_normal", 0.99
  )
)
```

## Measurements and a paired comparison

Squared error measures the error of each posterior mean, rather than the
mean squared distance of all posterior draws from the truth. Each
measurement row has a `target` and `value`. The runner adds study
identities, including `dataset_id` and `method_id`; the measurement name
distinguishes different quantities for the same target.

``` r

coefficient_error <- function(artifacts, context) {
  truth <- artifacts$data$truth
  means <- posterior::summarise_draws(artifacts$draws, "mean")
  estimate <- means$mean[match(names(truth), means$variable)]
  if (anyNA(estimate)) stop("A truth parameter has no posterior estimate")
  data.frame(target = names(truth), value = (estimate - truth)^2)
}

count_divergences <- function(artifacts, context) {
  counts <- artifacts$diagnostics$num_divergent
  if (is.null(counts)) stop("CmdStan returned no divergence counts")
  data.frame(target = "all_chains", value = sum(counts))
}

method_study <- with_measure(
  method_study, "squared_error", coefficient_error,
  needs = c("data", "draws")
)
method_study <- with_measure(
  method_study, "divergences", count_divergences, needs = "diagnostics"
)
```

The paired comparison consumes `test_scores` from two methods within one
`dataset_id`. It returns logit-normal minus Dirichlet total held-out log
score; a positive difference favors logit-normal for that dataset. It
does not combine training LOO estimates with held-out scores or pool
observations from different replicates.

``` r

compare_test_scores <- function(results, artifacts, context) {
  left <- artifacts[["logit_normal"]]$test_scores
  right <- artifacts[["dirichlet"]]$test_scores
  if (is.null(left) || is.null(right)) {
    return(data.frame(
      target = "logit_normal_minus_dirichlet",
      value = NA_real_, status = "missing_method_artifact"
    ))
  }
  if (!identical(left$observation_id, right$observation_id)) {
    stop("Methods scored different test observations")
  }
  if (!all(is.finite(c(left$log_score, right$log_score)))) {
    stop("A predictive score is not finite")
  }
  data.frame(
    target = "logit_normal_minus_dirichlet",
    value = sum(left$log_score - right$log_score), status = "complete"
  )
}

method_study <- with_comparison(
  method_study, "paired_test_log_score", compare_test_scores,
  by = "dataset_id", needs = "test_scores"
)
method_study <- with_retention(
  method_study, artifacts = c("draws", "test_scores", "diagnostics")
)
plan_study(method_study, replicates = 10L, seed = 2718L)
#> gdr2-prior-and-sampler
#> 2 conditions x 10 replicates = 20 datasets
#> 3 methods; 60 fits
#>                   name       needs                           stage
#>          squared_error data, draws                         per fit
#>            divergences diagnostics                         per fit
#>  paired_test_log_score test_scores within dataset, before eviction
#> Retain: draws, test_scores, diagnostics
#> Discard after use: data
#> Plan only: no data, models or cache records have been evaluated.
```

This declares 20 generated datasets and 60 method attempts. Attempt
records preserve failures; a missing method artifact produces an
explicit incomplete comparison. Divergent fits remain visible through
their diagnostics. No selection rule is applied here. The retained
scores support inspection of the paired calculation, and retained draws
support later parameter summaries.

``` r

method_run <- run_study(
  method_study, replicates = 10L, seed = 2718L,
  path = "runs/gdr2-method-example", workers = 1L
)
assessment <- assess_study(method_run)
```

## Further method variation

[GroupR2priors](https://github.com/jear2412/GroupR2priors/blob/bf8b361d97f6c49fe0be1fa242474fdb19d0307b/scripts/aux_functions/all_auxfunctions.R#L243-L315)
adds grouping assumptions and model-specific data preparation. Represent
those as method settings and fit functions, with named coefficient
truths retained in the generated data. Its parameter summaries also
distinguish null and non-null coefficients. That requires an explicit
target grouping in the measurement, beyond the per-parameter squared
errors shown here.

[Latent-Composite-HSGPs](https://github.com/jear2412/Latent-Composite-HSGPs/blob/c33a504c563720e4a19c3736272156df30666755/Simulation%20scripts/pcGP_data_scenario_script.R#L65-L160)
pregenerates data and compares exact GP and HSGP fits. Such methods
would supply their own rstan fit functions and chain-preserving
extractors, plus named latent truth vectors. Their approximation
settings belong to method identity. This vignette supplies neither those
fit functions nor the study’s calibration analysis; validating them
requires a separate case study with its own model sources and estimands.
