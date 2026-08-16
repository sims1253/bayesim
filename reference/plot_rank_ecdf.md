# Plot SBC rank ECDF with simultaneous confidence band

Plots the empirical CDF of SBC ranks against the uniform CDF (the
diagonal), with a simultaneous confidence band following Säilynoja,
Bürkner, and Vehtari (2022). The band is calibrated so that, under
correct calibration, the *entire* ECDF stays within it with probability
alpha; deviations anywhere along the band therefore indicate
miscalibration at level 1 - alpha.

Ranks are normalized per task: each task's ranks are scaled by that
task's own support, `(rank + 0.5) / n_ranks` with `n_ranks` = support +
1 (kept post-thinning draws + 1). When tasks in a panel have different
supports (e.g.
[`rank_metric()`](https://sims1253.github.io/bayesim/reference/RankMetric.md)
with `thin = "auto"` under autocorrelation), pooling on the panel
maximum would squash small-support tasks' ranks toward zero and
manufacture apparent miscalibration; per-task normalization avoids that
artifact, and a warning notes that the simultaneous band, which assumes
iid ranks on a common support, is then approximate. Legacy results
without a recorded `n_ranks` fall back per task to `n_draws`; with a
historical thinning stride \> 1 the true support is unknown, so that
fallback only bounds it.

## Usage

``` r
plot_rank_ecdf(ranks, alpha = 0.95, by = NULL)
```

## Arguments

- ranks:

  A tibble from
  [`sbc_ranks()`](https://sims1253.github.io/bayesim/reference/sbc_ranks.md),
  or a `bayesim_simulation_result`.

- alpha:

  Coverage level of the simultaneous confidence band (default 0.95).

- by:

  Optional character vector of condition columns to facet by. These
  columns are preserved by
  [`sbc_ranks()`](https://sims1253.github.io/bayesim/reference/sbc_ranks.md)
  for simulation results. Using `by` computes a separate ECDF and
  simultaneous band per condition cell instead of pooling ranks across
  cells.

## Value

A ggplot object.

## References

Säilynoja T, Bürkner PC, Vehtari A (2022). Graphical test for discrete
uniformity and its applications in goodness-of-fit evaluation.
*Statistics and Computing*, 32(2).

## Examples

``` r
if (FALSE) { # \dontrun{
plot_rank_ecdf(sbc_ranks(result))
plot_rank_ecdf(sbc_ranks(result), alpha = 0.99)
plot_rank_ecdf(sbc_ranks(result), by = "data_n")
} # }
```
