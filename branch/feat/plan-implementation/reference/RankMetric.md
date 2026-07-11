# SBC Rank Metric

Computes Simulation-Based Calibration ranks: for each parameter in
`vars_of_interest`, counts how many (possibly thinned) posterior draws
are below the data-generating `true_params` value. Under correct
calibration the ranks are uniformly distributed on `0..n_ranks-1`.

Autocorrelation in the posterior draws biases SBC rank uniformity (Talts
et al. 2018, §4.1). By default (`thin = "auto"`) the draws are thinned
toward the minimum bulk-ESS across the ranked variables before ranking:
this keeps ~ESS equally spaced draws, restoring near-independent samples
so the rank distribution is comparable to the standard SBC uniformity
test. `n_ranks` (posterior sample size after thinning + 1 possible
ranks) is reported per variable and is required by the SBC diagnostics
(`plot_rank_ecdf`).

Ranks use the strict comparison `draw < truth`, which is appropriate for
continuous posterior distributions where exact ties have probability
zero. For discrete parameters or parameters with boundary point masses,
use a custom metric with randomized tie-breaking.

Constructor for RankMetric.

## Usage

``` r
RankMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean",
  thin = "auto"
)

rank_metric(name = "rank", thin = "auto")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "rank".

- thin:

  Thinning policy: `"auto"` (default), `FALSE`, or an integer stride.

## Value

A `RankMetric` object.

A `RankMetric` object.

## Examples

``` r
rank_metric()
#> <bayesim::RankMetric>
#>  @ name        : chr "rank"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ thin        : chr "auto"
rank_metric(thin = FALSE)
#> <bayesim::RankMetric>
#>  @ name        : chr "rank"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ thin        : logi FALSE
```
