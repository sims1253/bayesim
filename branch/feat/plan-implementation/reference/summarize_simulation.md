# Aggregate simulation results per condition

Computes per-condition aggregates of the wide summary tibble: mean and
median of each metric column, Monte Carlo standard errors (MCSE),
replicate counts, and failure/convergence-failure rates. Returns a tidy
tibble with one row per condition.

Aggregation follows each metric's declared `summary_type` (E4; see
[Metric](https://sims1253.github.io/bayesim/reference/Metric.md)):
`"mean"` columns get a `sd / sqrt(n)` MCSE, `"proportion"` columns (e.g.
coverage) get `sqrt(p(1-p) / n)`, and `"none"` columns (e.g. SBC ranks)
are excluded from aggregation. Columns from unknown or user-defined
sources default to `"mean"`. MCSE formulas follow rsimsum (Gasparini,
2018).

## Usage

``` r
summarize_simulation(result, by = NULL, metrics = NULL)
```

## Arguments

- result:

  A `bayesim_simulation_result` object (uses `result$summary`), or a
  data.frame/tibble of per-task metrics. Passing the full result is
  preferred: it carries the metrics' `summary_type` declarations.

- by:

  Character vector of grouping columns (conditions). Defaults to the
  `data_*`/`fit_*` grid columns found in the summary plus any other
  non-numeric condition columns (excluding `task_id` and `status`).

- metrics:

  Character vector of metric columns to aggregate. Defaults to all
  numeric columns not in `by` and not metadata (`task_id`, `rep_idx`,
  `status`, `*timing*`).

## Value

A tibble with one row per condition: the `by` columns, then for each
metric `<m>_mean`, `<m>_median`, `<m>_sd`, `<m>_mcse`, plus `n_reps`,
`n_failed`, `failure_rate`.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_simulation(config, progress = FALSE)
summarize_simulation(result)
summarize_simulation(result, by = "model", metrics = "rmse__value")
} # }
```
