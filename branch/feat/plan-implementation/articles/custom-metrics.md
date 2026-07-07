# Custom Metrics

``` r

library(bayesim)
```

A metric computes one row’s worth of values per task — a summary of a
single fit against a single generated dataset. Aggregation *across*
tasks (bias, coverage, MCSEs) happens later, in
[`summarize_simulation()`](https://sims1253.github.io/bayesim/reference/summarize_simulation.md)
and
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md).
This vignette covers writing your own metric: the contract, the output
schema, `summary_type`, and how large outputs are externalized.

## The Metric contract

Extend the `Metric` S7 class and implement
[`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md):

``` r

MADMetric <- S7::new_class(
  "MADMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "mad"),
    needs = S7::new_property(S7::class_character, default = "predictions"),
    required = S7::new_property(S7::class_logical, default = FALSE)
  )
)

S7::method(compute_metric, MADMetric) <- function(
  metric, fit_result, data_bundle, context, task_ctx
) {
  if (is.null(context$predictions) || is.null(data_bundle$test)) {
    return(list(value = NA_real_))
  }
  actual <- data_bundle$test[[data_bundle$response]]
  predicted <- context$predictions$predicted_mean
  list(value = stats::median(abs(predicted - actual)))
}

# A constructor keeps configs readable.
mad_metric <- function(name = "mad") MADMetric(name = name)
```

[`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md)
receives:

- `fit_result` — the `bayesim_fit_result` (with `$draws`,
  `$diagnostics`, and `$fit` if retained).
- `data_bundle` — the generator’s output (`$train`, `$test`,
  `$response`, `$true_params`, `$vars_of_interest`).
- `context` — precomputed values, driven by your metric’s `needs`
  property: `"predictions"` (a
  [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
  result), `"log_lik"` (an S x N matrix), `"loo"` (elpd summary plus
  PSIS machinery). Declare only what you use; the engine computes each
  at most once per task and shares it across metrics. Predictions and
  log-lik are evaluated on the **test set** when one exists, otherwise
  on the training data.
- `task_ctx` — task identity (`task_id`, `data_idx`, `fit_idx`,
  `rep_idx`, `seed`).

Degrade to `NA` when inputs are missing (as above) rather than erroring:
metrics with `required = FALSE` (the default) record NA on failure,
while `required = TRUE` metrics fail the whole task.

## The output schema

[`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md)
must return a named list where every element is either a scalar atomic
value or a **named** numeric vector. No matrices, data frames, or nested
lists — task results must stay flat and cheap to store.

The engine flattens the output into summary columns as
`<metric_name>__<field>`, and named vectors expand per element:

``` r

list(value = 0.3, n_obs = 50L)
#> mad__value, mad__n_obs

list(by_param = c(Intercept = 1, x = 0))
#> mad__by_param__Intercept, mad__by_param__x
```

[`validate_metric_output()`](https://sims1253.github.io/bayesim/reference/validate_metric_output.md)
enforces this schema; run your metric once against a fixture fit in your
tests and pass the result through it.

## summary_type: how your metric aggregates

[`summarize_simulation()`](https://sims1253.github.io/bayesim/reference/summarize_simulation.md)
needs to know how to aggregate each metric’s columns across replicates.
Declare it with the `summary_type` property:

- `"mean"` (default) — report mean/median/sd with an `sd/sqrt(n)` MCSE.
- `"proportion"` — for 0/1 outcomes like coverage; MCSE is
  `sqrt(p(1-p)/n)`.
- `"none"` — never aggregate (e.g. SBC ranks, which are analyzed as a
  distribution via
  [`sbc_ranks()`](https://sims1253.github.io/bayesim/reference/sbc_ranks.md),
  not averaged).

``` r

HitMetric <- S7::new_class(
  "HitMetric",
  parent = Metric,
  properties = list(
    name = S7::new_property(S7::class_character, default = "hit"),
    summary_type = S7::new_property(S7::class_character, default = "proportion")
  )
)
```

## Parallel safety

Under `run_simulation(config, workers = N)`, metric objects are crated
and shipped to daemon processes. Keep
[`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md)
methods self-contained: call package functions by namespace
([`stats::median`](https://rdrr.io/r/stats/median.html), not a
re-exported alias) and do not reference variables captured from your
interactive session.

## Externalization of large outputs

Named numeric vectors longer than 50 elements (or above 64 kB) are not
inlined into the summary — with a `result_path` set, they are written to
`<result_path>/artifacts/metrics/` and the summary records a pointer
(`<metric>__<field>__artifact_path`, `__artifact_hash`, `__n_values`).
This keeps summaries navigable when a metric emits, say, per-observation
Pareto-k values across ten thousand tasks. Without a `result_path`,
large vectors are inlined as usual.

## Using a custom metric

``` r

gen <- function(data_spec, task_ctx) {
  n <- data_spec$n
  x <- stats::rnorm(2 * n)
  y <- 2 * x + stats::rnorm(2 * n)
  idx <- seq_len(n)
  list(
    train = data.frame(y = y, x = x)[idx, ],
    test = data.frame(y = y, x = x)[-idx, ],
    response = "y",
    true_params = c(Intercept = 0, x = 2, sigma = 1),
    vars_of_interest = c("Intercept", "x", "sigma")
  )
}

config <- simulation_config(
  data_grid = data.frame(n = 50),
  fit_grid = data.frame(model = "linear"),
  data_generator = gen,
  fitter = LinearRegressionFitter(n_draws = 200L),
  metrics = list(mad_metric(), posterior_summary_metric()),
  n_replicates = 4L,
  seed = 42L
)

result <- run_simulation(config, progress = FALSE)
#> 4 tasks = 1 data x 1 fit x 4 reps
#> ℹ Starting simulation with 4 tasks
summarize_simulation(result, metrics = "mad__value")
#> # A tibble: 1 × 9
#>   data_n fit_model n_reps n_failed failure_rate mad__value_mean
#>    <dbl> <chr>      <int>    <int>        <dbl>           <dbl>
#> 1     50 linear         4        0            0              NA
#> # ℹ 3 more variables: mad__value_median <dbl>, mad__value_sd <dbl>,
#> #   mad__value_mcse <dbl>
```

## Next steps

- [`vignette("custom-fitters")`](https://sims1253.github.io/bayesim/articles/custom-fitters.md)
  — the other half of the extension surface.
- [`?Metric`](https://sims1253.github.io/bayesim/reference/Metric.md) —
  the canonical contract reference.
- [`vignette("design-of-simulation-studies")`](https://sims1253.github.io/bayesim/articles/design-of-simulation-studies.md)
  — which performance measures to compute across tasks, and their MCSEs.
