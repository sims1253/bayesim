# Experimental grammar: declare, run, reuse

Describe the experiment with ordinary R functions. The runner supplies
shared data, deterministic random streams, saved results and comparison
groups. This experimental interface is sourced explicitly; the supported
package API still starts with
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md).

``` text
condition × replicate → one dataset → each method
                                       ↓
                              extract → measure → compare
                                       ↓             ↓
                              retained artifacts + results
                                       ↓
                              new measures and assessment
```

This example runs a small normal-mean study. Observations have known
unit variance. Two methods differ only in their normal prior’s scale.
Their posteriors are analytic, so this vignette needs no sampler.

``` r

normal_method <- function(prior_sd) {
  study_method(
    fit = function(data, context) {
      precision <- length(data) + 1 / context$settings$prior_sd^2
      list(mean = sum(data) / precision, sd = sqrt(1 / precision))
    },
    extract = list(draws = function(fit, data, context) {
      rnorm(400, fit$mean, fit$sd)
    }),
    settings = list(prior_sd = prior_sd)
  )
}

experiment <- study(
  "normal-mean",
  conditions = data.frame(condition_id = c("small", "large"), n = c(10L, 50L)),
  generate = function(condition, context) rnorm(condition$n, mean = 1)
) |>
  with_methods(regularized = normal_method(1), diffuse = normal_method(10)) |>
  with_measure("estimate", function(artifacts, context) {
    data.frame(target = "mean", value = mean(artifacts$draws))
  }, needs = "draws") |>
  with_retention("draws")

plan_study(experiment, replicates = 3)
#> normal-mean
#> 2 conditions x 3 replicates = 6 datasets
#> 2 methods; 12 fits
#>      name needs   stage
#>  estimate draws per fit
#> Retain: draws
#> Discard after use: 
#> Plan only: no data, models or cache records have been evaluated.
```

The draws here are independent posterior samples. An MCMC extractor
should preserve chains and retain any diagnostics needed for later
selection.

``` r

run_path <- tempfile("normal-mean-")
first <- run_study(experiment, replicates = 3, path = run_path)
first$measurements[, c("condition_id", "replicate", "method_id", "value")]
#>    condition_id replicate   method_id     value
#> 1         small         1 regularized 1.1332277
#> 2         small         1     diffuse 1.2724442
#> 3         large         1 regularized 0.9444568
#> 4         large         1     diffuse 0.9634327
#> 5         small         2 regularized 0.6258363
#> 6         small         2     diffuse 0.6751752
#> 7         large         2 regularized 0.8710943
#> 8         large         2     diffuse 0.8925380
#> 9         small         3 regularized 0.6224212
#> 10        small         3     diffuse 0.6804151
#> 11        large         3 regularized 1.0074941
#> 12        large         3     diffuse 1.0476324
```

Add a 90% interval and two replicates. The original fits stay cached;
the new measurement uses their saved draws. New replicates get new
datasets and fits.

``` r

extended <- experiment |>
  with_measure("interval90", function(artifacts, context) {
    interval <- quantile(artifacts$draws, c(0.05, 0.95))
    data.frame(target = "mean", lower = interval[[1]], upper = interval[[2]])
  }, needs = "draws")

later <- run_study(extended, replicates = 5, path = run_path)
stopifnot(identical(first$attempts, later$attempts[later$attempts$replicate <= 3, ]))
table(later$measurements$measurement)
#> 
#>   estimate interval90 
#>         20         20
```

This reuse requires unchanged scientific dependencies. Cache keys track
function bodies, settings and explicit versions, but cannot detect
changed global helpers, package versions or external files. Increment
the relevant declaration’s `version` when those dependencies change.

The case studies test the boundaries of this interface:

| Study | What the declaration must express |
|----|----|
| [Choice of likelihood](https://sims1253.github.io/bayesim/articles/experimental-likelihood-study.md) | Full candidate grids, shared datasets, predictive comparisons and changing eligibility |
| [Primed priors](https://sims1253.github.io/bayesim/articles/experimental-primed-priors.md) | Preparation shared across inner replicates and data reused in refits |
| [Method comparison](https://sims1253.github.io/bayesim/articles/experimental-method-comparison.md) | Direct CmdStan methods, sampler settings and pointwise comparison inputs |

Those vignettes build plans without running their expensive fits. This
runner is a reference implementation: individual RDS records and
in-memory result tables still need replacement or validation for a study
with millions of fits. It also requires exclusive access to a run
directory. Its purpose is to test whether the study declaration captures
the scientific decisions before committing to a production execution and
storage design.
