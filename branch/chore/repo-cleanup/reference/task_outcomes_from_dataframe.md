# Migrate the legacy flat checkpoint view to canonical task outcomes.

The flat view produced by
[`results_to_dataframe()`](https://sims1253.github.io/bayesim/reference/results_to_dataframe.md)
lays columns out in a fixed order: `task_id`, `status`, `stop_reason`,
metric columns, `truth__*` columns, diagnostic columns, then the
optional `error_class`/`error_message`/`timing_total` trailer. Metric
and diagnostic columns share a namespace, so the two groups are
separated positionally around the `truth__` block. When a row carries no
truth block the groups are contiguous and indistinguishable; every value
column is then restored as a metric so the flat representation still
round-trips unchanged.

## Usage

``` r
task_outcomes_from_dataframe(results_df)
```
