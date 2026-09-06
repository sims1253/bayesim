# Getting Started with bayesim

This vignette defines a data generator, runs a Bayesian linear
regression study, and summarizes its results. It uses
[`LinearRegressionFitter()`](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
for exact conjugate inference without Stan. For Stan models, use
[`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
or
[`CmdStanFitter()`](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md).

## Run a study

``` r

library(bayesim)
```

### Generate data

A data generator is a function with the signature
`(data_spec, task_ctx)` returning a named list. `data_spec` contains one
row of `data_grid`; `task_ctx` identifies the task. Use R’s random
functions normally. Set the study seed in
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md),
not inside the generator.

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

### Choose conditions and metrics

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

### Run the simulation

``` r

result <- run_simulation(config, progress = FALSE)
#> 8 tasks = 2 data x 1 fit x 4 reps
#> ℹ Starting simulation with 8 tasks
#> 
#> ✔ Simulation complete: 8/8 tasks succeeded in 0.2s
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
#>   Total time: 0.17 s
```

### Inspect task results

``` r

result$summary[c("task_id", "data_n", "status", "posterior_summary__mean__x")]
#>            task_id data_n  status posterior_summary__mean__x
#> 1 d001_f001_r00001     50 success                   1.854738
#> 2 d001_f001_r00002     50 success                   1.859731
#> 3 d001_f001_r00003     50 success                   2.064685
#> 4 d001_f001_r00004     50 success                   1.760973
#> 5 d002_f001_r00001    100 success                   2.035980
#> 6 d002_f001_r00002    100 success                   1.842911
#> 7 d002_f001_r00003    100 success                   2.177546
#> 8 d002_f001_r00004    100 success                   1.881039
```

Each row is one task. Columns include `task_id`, `status`,
`timing_total`, the grid columns, and one column per metric field
(`posterior_summary__mean__x`, `truth__x`, …). The `truth__*` columns
record the data-generating truth for parameter-recovery analysis.

## Summarize results

Each sample size is a separate condition. Group by `data_n` to keep
their results separate. Select the slope’s posterior mean to keep the
table small:

``` r

summarize_simulation(
  result,
  by = "data_n",
  metrics = "posterior_summary__mean__x"
)
#> # A tibble: 2 × 9
#>   data_n n_reps n_failed failure_rate posterior_summary__mean__x_n_used
#>    <dbl>  <int>    <int>        <dbl>                             <int>
#> 1     50      4        0            0                                 4
#> 2    100      4        0            0                                 4
#> # ℹ 4 more variables: posterior_summary__mean__x_mean <dbl>,
#> #   posterior_summary__mean__x_median <dbl>,
#> #   posterior_summary__mean__x_sd <dbl>, posterior_summary__mean__x_mcse <dbl>
```

The output reports the mean, median, standard deviation, and Monte Carlo
standard error (MCSE) across replicates, plus task and failure counts.
Use `metric_cols(result, "posterior_summary")` to find other metric
columns.

## Performance measures

Use
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)
for bias, empirical SE, MSE, coverage, and average model SE, each with
its MCSE (Morris, White & Crowther, 2019):

``` r

pm <- performance_measures(result, estimand = "x", by = "data_n")
pm
#> # A tibble: 12 × 7
#>    data_n estimand measure    value     mcse n_sim truth_mode
#>     <dbl> <chr>    <chr>      <dbl>    <dbl> <int> <chr>     
#>  1     50 x        bias     -0.115   0.0640      4 fixed     
#>  2     50 x        emp_se    0.128   0.0523      4 fixed     
#>  3     50 x        mse       0.0255  0.0112      4 fixed     
#>  4     50 x        model_se  0.140   0.00704     4 fixed     
#>  5     50 x        coverage  1       0           4 fixed     
#>  6     50 x        n_sim     4      NA           4 fixed     
#>  7    100 x        bias     -0.0156  0.0767      4 fixed     
#>  8    100 x        emp_se    0.153   0.0627      4 fixed     
#>  9    100 x        mse       0.0179  0.00659     4 fixed     
#> 10    100 x        model_se  0.0996  0.00489     4 fixed     
#> 11    100 x        coverage  1       0           4 fixed     
#> 12    100 x        n_sim     4      NA           4 fixed
```

This pairs `truth__x` with the per-task posterior summary to give, per
condition cell, the bias and calibration of the estimator.

## Run a larger study

This example uses four replicates per condition to keep it quick. Choose
more replicates for a precise estimate; the [study design
guide](https://sims1253.github.io/bayesim/articles/design-of-simulation-studies.md)
explains how MCSE helps you decide.

Add `workers = 4` to
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
to use four parallel workers. For long runs, set `result_path` and
`checkpoint_every` when constructing the configuration. Resume
interrupted work with
`resume_simulation("my_simulation", config = config)`, using the
original configuration. Script-defined generators require that
configuration because the saved manifest cannot reconstruct their
closures.

See [parallel runs and
checkpoints](https://sims1253.github.io/bayesim/articles/parallel-and-hpc.md)
for a complete example, and
[reproducibility](https://sims1253.github.io/bayesim/articles/reproducibility.md)
for seed and resume requirements.

## Next steps

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
