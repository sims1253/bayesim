# Getting Started with bayesim

## Introduction

bayesim is a simulation framework for Bayesian modeling studies. It
provides:

- **Reproducible execution** across sequential, parallel (mirai), and
  resumed runs
- **A model bank** that compiles each distinct Stan model once and
  reuses it across tasks, eliminating per-task recompilation
- **Checkpoint/resume** capabilities for long-running simulations
- **Memory-bounded execution** with configurable artifact retention
- **Typed generators, metrics, and an analysis layer** – no hand-rolled
  dplyr

This vignette runs end-to-end with
[`LinearRegressionFitter()`](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
– exact conjugate Bayesian linear regression (real posteriors, no Stan,
milliseconds per fit). For Stan/brms models, swap it for
[`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
or
[`CmdStanFitter()`](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md).

## Basic Usage

### Setup

``` r

library(bayesim)
```

### Define a Data Generator

A data generator is a function with the signature
`(data_spec, task_ctx)` returning a named list. It must consume the
**ambient** RNG state (the engine restores a per-task L’Ecuyer stream
before each call), so do not re-seed inside.

``` r

my_data_generator <- function(data_spec, task_ctx) {
  n <- data_spec$n
  x <- stats::rnorm(n)
  y <- data_spec$intercept + data_spec$slope * x +
    stats::rnorm(n, sd = data_spec$sigma)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(
      Intercept = data_spec$intercept,
      x = data_spec$slope,
      sigma = data_spec$sigma
    ),
    vars_of_interest = c("Intercept", "x", "sigma")
  )
}
```

For common patterns (fixed truth, prior-predictive, inverse forward
sampling) use the factory constructors
[`fixed_truth_generator()`](https://sims1253.github.io/bayesim/reference/fixed_truth_generator.md),
[`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md),
[`ifs_generator()`](https://sims1253.github.io/bayesim/reference/ifs_generator.md).

### Create a Configuration

``` r

config <- simulation_config(
  data_grid = data.frame(
    n = c(50, 100),
    intercept = 1,
    slope = 2,
    sigma = 1
  ),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_data_generator,
  fitter = LinearRegressionFitter(n_draws = 500L),
  metrics = list(
    posterior_summary_metric(),
    sampler_diagnostics_metric()
  ),
  n_replicates = 4L,
  seed = 42L
)
```

### Run the Simulation

``` r

result <- run_simulation(config, progress = FALSE)
#> 8 tasks = 2 data x 1 fit x 4 reps
#> ℹ Starting simulation with 8 tasks
print(result)
#> <bayesim_simulation_result>
#>   Config fingerprint: 0443023694ffb037043507f14667a212f03a59b2d3913eb62fdbba651afdb417 
#>   Tasks: 8 
#>     - Success: 8 
#>     - Failed: 0 
#>     - Pending: 0 
#>     - Skipped (policy-stopped): 0 
#>   Metrics: posterior_summary__mean__Intercept, posterior_summary__mean__x, posterior_summary__mean__sigma, posterior_summary__median__Intercept, posterior_summary__median__x, posterior_summary__median__sigma  ... 
#>   Task grid: 8 rows x 7 cols
#>   Total time: 0.22 s
```

### Examine Results

``` r

head(result$summary)
#>            task_id  status stop_reason posterior_summary__mean__Intercept
#> 1 d001_f001_r00001 success        <NA>                          0.8448899
#> 2 d001_f001_r00002 success        <NA>                          1.0002337
#> 3 d001_f001_r00003 success        <NA>                          0.9888108
#> 4 d001_f001_r00004 success        <NA>                          0.9185220
#> 5 d002_f001_r00001 success        <NA>                          0.9700026
#> 6 d002_f001_r00002 success        <NA>                          1.1828572
#>   posterior_summary__mean__x posterior_summary__mean__sigma
#> 1                   1.854738                      0.8546019
#> 2                   1.859731                      0.9285260
#> 3                   2.064685                      1.0919644
#> 4                   1.760973                      0.9530450
#> 5                   2.035980                      0.9592429
#> 6                   1.842911                      1.0116062
#>   posterior_summary__median__Intercept posterior_summary__median__x
#> 1                            0.8479619                     1.844666
#> 2                            1.0013767                     1.858353
#> 3                            0.9896985                     2.068306
#> 4                            0.9151737                     1.765062
#> 5                            0.9671287                     2.040978
#> 6                            1.1837680                     1.848306
#>   posterior_summary__median__sigma posterior_summary__sd__Intercept
#> 1                        0.8478846                       0.13191406
#> 2                        0.9217226                       0.14238030
#> 3                        1.0888450                       0.15523026
#> 4                        0.9467269                       0.14825666
#> 5                        0.9598104                       0.09811934
#> 6                        1.0003037                       0.10058448
#>   posterior_summary__sd__x posterior_summary__sd__sigma
#> 1                0.1291764                   0.08570589
#> 2                0.1266755                   0.08849451
#> 3                0.1554954                   0.11115861
#> 4                0.1477083                   0.08892455
#> 5                0.1016660                   0.06865704
#> 6                0.1125469                   0.07563001
#>   posterior_summary__q_lower__Intercept posterior_summary__q_lower__x
#> 1                             0.5958407                      1.613671
#> 2                             0.7196671                      1.620410
#> 3                             0.6995356                      1.765564
#> 4                             0.6204440                      1.485925
#> 5                             0.7918690                      1.839979
#> 6                             0.9750488                      1.607716
#>   posterior_summary__q_lower__sigma posterior_summary__q_upper__Intercept
#> 1                         0.7197477                              1.086398
#> 2                         0.7673389                              1.277838
#> 3                         0.8988157                              1.289214
#> 4                         0.8099143                              1.196901
#> 5                         0.8309318                              1.161128
#> 6                         0.8914731                              1.378350
#>   posterior_summary__q_upper__x posterior_summary__q_upper__sigma
#> 1                      2.125122                          1.050315
#> 2                      2.115247                          1.106836
#> 3                      2.382265                          1.316516
#> 4                      2.042422                          1.148037
#> 5                      2.230406                          1.104865
#> 6                      2.060749                          1.168051
#>   sampler_diagnostics__rhat_max sampler_diagnostics__ess_bulk_min
#> 1                             1                                NA
#> 2                             1                                NA
#> 3                             1                                NA
#> 4                             1                                NA
#> 5                             1                                NA
#> 6                             1                                NA
#>   sampler_diagnostics__ess_tail_min sampler_diagnostics__divergent
#> 1                                NA                              0
#> 2                                NA                              0
#> 3                                NA                              0
#> 4                                NA                              0
#> 5                                NA                              0
#> 6                                NA                              0
#>   sampler_diagnostics__max_treedepth truth__Intercept truth__x truth__sigma
#> 1                                  0                1        2            1
#> 2                                  0                1        2            1
#> 3                                  0                1        2            1
#> 4                                  0                1        2            1
#> 5                                  0                1        2            1
#> 6                                  0                1        2            1
#>   rhat_max ess_bulk ess_tail divergent max_treedepth timing_total rep_idx
#> 1        1      500      500         0             0  0.015394926       1
#> 2        1      500      500         0             0  0.034210920       2
#> 3        1      500      500         0             0  0.037766695       3
#> 4        1      500      500         0             0  0.004510880       4
#> 5        1      500      500         0             0  0.004056454       1
#> 6        1      500      500         0             0  0.004071951       2
#>   data_n data_intercept data_slope data_sigma fit_model
#> 1     50              1          2          1    linear
#> 2     50              1          2          1    linear
#> 3     50              1          2          1    linear
#> 4     50              1          2          1    linear
#> 5    100              1          2          1    linear
#> 6    100              1          2          1    linear
```

Each row is one task. Columns include `task_id`, `status`,
`timing_total`, the grid columns, and one column per metric field
(`posterior_summary__mean__x`, `truth__x`, …). The `truth__*` columns
record the data-generating truth (E1), enabling parameter-recovery
analysis.

## Summarizing Results

[`summarize_simulation()`](https://sims1253.github.io/bayesim/reference/summarize_simulation.md)
aggregates the per-task summary by condition, with Monte Carlo standard
errors (MCSE):

``` r

agg <- summarize_simulation(result, by = "data_n")
agg
#> # A tibble: 2 × 154
#>   data_n n_reps n_failed failure_rate posterior_summary__mean__Intercept_n_used
#>    <dbl>  <int>    <int>        <dbl>                                     <int>
#> 1     50      4        0            0                                         4
#> 2    100      4        0            0                                         4
#> # ℹ 149 more variables: posterior_summary__mean__Intercept_mean <dbl>,
#> #   posterior_summary__mean__Intercept_median <dbl>,
#> #   posterior_summary__mean__Intercept_sd <dbl>,
#> #   posterior_summary__mean__Intercept_mcse <dbl>,
#> #   posterior_summary__mean__x_n_used <int>,
#> #   posterior_summary__mean__x_mean <dbl>,
#> #   posterior_summary__mean__x_median <dbl>, …
```

For each metric you get `<metric>_mean`, `<metric>_median`,
`<metric>_sd`, and `<metric>_mcse`, plus `n_reps`, `n_failed`, and
`failure_rate`.

Studies with several metrics can produce wide summaries – each metric
contributes `_n_used`, `_mean`, `_median`, `_sd`, and `_mcse` columns,
so the aggregate can easily exceed 100 columns. Nothing is dropped or
truncated; narrow the aggregation with the `metrics` argument, and use
[`metric_cols()`](https://sims1253.github.io/bayesim/reference/metric_cols.md)
to list the flattened columns of a single metric:

``` r

metric_cols(result, "posterior_summary")
#>                         mean__Intercept                                 mean__x 
#>    "posterior_summary__mean__Intercept"            "posterior_summary__mean__x" 
#>                             mean__sigma                       median__Intercept 
#>        "posterior_summary__mean__sigma"  "posterior_summary__median__Intercept" 
#>                               median__x                           median__sigma 
#>          "posterior_summary__median__x"      "posterior_summary__median__sigma" 
#>                           sd__Intercept                                   sd__x 
#>      "posterior_summary__sd__Intercept"              "posterior_summary__sd__x" 
#>                               sd__sigma                      q_lower__Intercept 
#>          "posterior_summary__sd__sigma" "posterior_summary__q_lower__Intercept" 
#>                              q_lower__x                          q_lower__sigma 
#>         "posterior_summary__q_lower__x"     "posterior_summary__q_lower__sigma" 
#>                      q_upper__Intercept                              q_upper__x 
#> "posterior_summary__q_upper__Intercept"         "posterior_summary__q_upper__x" 
#>                          q_upper__sigma 
#>     "posterior_summary__q_upper__sigma"
agg_post <- summarize_simulation(
  result,
  by = "data_n",
  metrics = metric_cols(result, "posterior_summary", fields = "mean")
)
```

## Performance measures

For estimator-performance measures in the sense of Morris, White &
Crowther (2019) – bias, empirical SE, MSE, coverage, average model SE,
each with its MCSE – use
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md):

``` r

pm <- performance_measures(result, estimand = "x")
pm
#> # A tibble: 12 × 11
#>    data_n data_intercept data_slope data_sigma fit_model estimand measure 
#>     <dbl>          <dbl>      <dbl>      <dbl> <chr>     <chr>    <chr>   
#>  1     50              1          2          1 linear    x        bias    
#>  2     50              1          2          1 linear    x        emp_se  
#>  3     50              1          2          1 linear    x        mse     
#>  4     50              1          2          1 linear    x        model_se
#>  5     50              1          2          1 linear    x        coverage
#>  6     50              1          2          1 linear    x        n_sim   
#>  7    100              1          2          1 linear    x        bias    
#>  8    100              1          2          1 linear    x        emp_se  
#>  9    100              1          2          1 linear    x        mse     
#> 10    100              1          2          1 linear    x        model_se
#> 11    100              1          2          1 linear    x        coverage
#> 12    100              1          2          1 linear    x        n_sim   
#> # ℹ 4 more variables: value <dbl>, mcse <dbl>, n_sim <int>, truth_mode <chr>
```

This pairs `truth__x` with the per-task posterior summary to give, per
condition cell, the bias and calibration of the estimator. It is the
function a methods paper actually needs.

## Parallel Execution (mirai)

Parallelism is owned by the engine and uses
[mirai](https://mirai.r-lib.org/) via
[purrr](https://purrr.tidyverse.org/)’s `in_parallel()` integration. The
simple path is the `workers` argument:

``` r

result <- run_simulation(config, progress = FALSE, workers = 4)
# daemons are set up for the run and torn down on exit
```

For the advanced/HPC path (remote daemons, TLS), call
[`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
yourself before
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md).
`SimulationConfig` accepts a `daemon_setup` function run once per daemon
(e.g. to configure the cmdstan path when using `BrmsFitter` on a
cluster).

## Checkpointing and Resume

For long-running simulations, write checkpoints to `result_path`:

``` r

config <- simulation_config(
  data_grid = data.frame(
    n = c(50, 100),
    intercept = 1,
    slope = 2,
    sigma = 1
  ),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_data_generator,
  fitter = LinearRegressionFitter(n_draws = 200L),
  metrics = list(posterior_summary_metric()),
  n_replicates = 1000L,
  seed = 42L,
  result_path = "my_simulation",
  checkpoint_every = 50L
)
result <- run_simulation(config, resume = "auto")
```

If interrupted, resume with the original config:

``` r

result <- resume_simulation("my_simulation", config = config)
```

Passing the config is required here because `my_data_generator` is a
closure defined in your script; the run manifest cannot rehydrate
script-defined closures. The configless form
`resume_simulation("my_simulation")` works only when every generator,
fitter, and metric is a namespaced package function or class.

## Metrics

The built-in metric library covers the standard surface:
`pred_bias_metric`, `pred_rmse_metric`, `pred_mae_metric`,
`pred_mse_metric`, `coverage_metric`, `pos_prob_metric`,
`posterior_summary_metric`, `sampler_diagnostics_metric` (convergence
and sampler diagnostics), `rank_metric` (SBC ranks), and LOO/test-set
variants (`elpd_loo_metric`, `rmse_loo_metric`, `r2_loo_metric`,
`elpd_test_metric`, `r2_test_metric`).

To write your own metric, extend the `Metric` class and implement
[`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md);
[`vignette("custom-metrics")`](https://sims1253.github.io/bayesim/articles/custom-metrics.md)
covers the full contract, the output schema, and how large outputs are
externalized.

## Reproducibility

Same seed + same package/backend versions + same platform produces
identical results. Across platforms (or Stan version changes) results
are statistically equivalent but may differ in the least significant
bits due to floating-point and backend differences.

## Next Steps

- [`vignette("design-of-simulation-studies")`](https://sims1253.github.io/bayesim/articles/design-of-simulation-studies.md)
  for aims, estimands, and choosing the number of replicates
- [`vignette("sbc-and-calibration")`](https://sims1253.github.io/bayesim/articles/sbc-and-calibration.md)
  for the calibration workflow
- [`vignette("brms-studies")`](https://sims1253.github.io/bayesim/articles/brms-studies.md)
  for Stan-backed studies and
  [`model_grid()`](https://sims1253.github.io/bayesim/reference/model_grid.md)
- [`vignette("custom-fitters")`](https://sims1253.github.io/bayesim/articles/custom-fitters.md)
  /
  [`vignette("custom-metrics")`](https://sims1253.github.io/bayesim/articles/custom-metrics.md)
  for extending the framework
- [`vignette("parallel-and-hpc")`](https://sims1253.github.io/bayesim/articles/parallel-and-hpc.md)
  and
  [`vignette("reproducibility")`](https://sims1253.github.io/bayesim/articles/reproducibility.md)
  for large runs and the determinism guarantees
