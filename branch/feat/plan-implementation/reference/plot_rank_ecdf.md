# Plot SBC rank ECDF with simultaneous confidence band

Plots the empirical CDF of SBC ranks against the uniform CDF (the
diagonal), with a simultaneous confidence band following Säilynoja,
Bürkner, and Vehtari (2022). The band is calibrated so that, under
correct calibration, the *entire* ECDF stays within it with probability
alpha; deviations anywhere along the band therefore indicate
miscalibration at level 1 - alpha.

## Usage

``` r
plot_rank_ecdf(ranks, alpha = 0.95)
```

## Arguments

- ranks:

  A tibble from
  [`sbc_ranks()`](https://sims1253.github.io/bayesim/reference/sbc_ranks.md),
  or a `bayesim_simulation_result`.

- alpha:

  Coverage level of the simultaneous confidence band (default 0.95).

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
} # }
```
