# Extract SBC ranks from a simulation result

Collects the per-task `rank__by_param` entries (from
[`rank_metric()`](https://sims1253.github.io/bayesim/reference/RankMetric.md))
into a long tibble with columns `task_id`, `param`, `rank`, `n_draws`.
Returns an empty tibble if no rank metric was computed.

## Usage

``` r
sbc_ranks(result)
```

## Arguments

- result:

  A `bayesim_simulation_result`.

## Value

A tibble with columns `task_id`, `param`, `rank`, `n_draws`.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- run_simulation(config, progress = FALSE)
ranks <- sbc_ranks(result)
} # }
```
