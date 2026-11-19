# Plot SBC rank histograms

One histogram panel per parameter. Requires ggplot2.

## Usage

``` r
plot_rank_hist(ranks)
```

## Arguments

- ranks:

  A tibble from
  [`sbc_ranks()`](https://sims1253.github.io/bayesim/reference/sbc_ranks.md),
  or a `bayesim_simulation_result`.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_rank_hist(sbc_ranks(result))
} # }
```
