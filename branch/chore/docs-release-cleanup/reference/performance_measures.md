# Simulation-method performance measures with Monte-Carlo standard errors

For condition cells with fixed data-generating truth, computes the
Morris, White & Crowther (2019, *Stat Med*) estimator-performance
measures: **bias**, **empirical SE**, **MSE**, **coverage**, **average
model SE**, and `n_sim`, each with its MCSE. When truth varies across
replicates (as in a prior-predictive SBC study), the fixed-truth Morris
names are not returned. Instead, the error distribution is described by
`mean_error`, `error_sd`, and `error_mse`; coverage and average model SE
remain available.

For each parameter the function pairs the data-generating
`truth__<param>` column with the per-task
`posterior_summary__*__<param>` columns (point estimate `mean`,
posterior `sd`, and interval `q_lower`/`q_upper`). Coverage uses the
interval when available; otherwise it falls back to a
`coverage__by_param__<param>` column if present.

Fixed-truth MCSE formulas follow Morris et al. / rsimsum:

- bias MCSE = sd(est - truth) / sqrt(n)

- empSE MCSE = sd / sqrt(2(n-1))

- MSE MCSE = sqrt(Var((est-truth)^2) / n)

- coverage MCSE = sqrt(p(1-p) / n)

- modelSE MCSE = sd(posterior_sd) / sqrt(n)

The bias MCSE uses the sd of the estimation errors `est - truth`, which
is valid under fixed and varying truth alike (with fixed truth the truth
is constant, so sd(est - truth) = sd(est)). For varying truth,
`mean_error` uses `sd(est-truth) / sqrt(n)`, `error_sd` uses
`error_sd / sqrt(2(n-1))`, and `error_mse` uses
`sqrt(Var((est-truth)^2) / n)`.

## Usage

``` r
performance_measures(
  result,
  estimand = NULL,
  estimator = c("mean", "median"),
  by = NULL
)
```

## Arguments

- result:

  A `bayesim_simulation_result` (uses `$summary`), or a data.frame of
  per-task metrics.

- estimand:

  Optional character; a single parameter name. When NULL (default), all
  parameters with both a `truth__*` and `posterior_summary__mean__*`
  column are analyzed.

- estimator:

  Character scalar naming the per-task point estimate to use: one of
  `"mean"` (default), `"median"`. Selects the corresponding
  `posterior_summary__<estimator>__<param>` column.

- by:

  Character vector of condition/grouping columns. Defaults to the
  data_grid and fit_grid columns found in the summary (excluding
  task_id, rep_idx, status, and metric columns).

## Value

A tidy tibble with columns: the `by` columns, `estimand`, `measure`,
`value`, `mcse`, `n_sim`, `truth_mode`. One row per condition x estimand
x measure. Fixed-truth measures use `bias`, `emp_se`, and `mse`;
varying-truth measures use `mean_error`, `error_sd`, and `error_mse`.
Both modes may also include `coverage`, `model_se`, and `n_sim`.

## References

Morris, White & Crowther (2019), *Using simulation studies to evaluate
statistical methods*, Statistics in Medicine.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_simulation(config, progress = FALSE)
performance_measures(result)
performance_measures(result, estimand = "x", by = "model")
} # }
```
