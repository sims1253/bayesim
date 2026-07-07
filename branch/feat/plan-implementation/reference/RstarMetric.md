# R\* Convergence Metric

Computes the R\* MCMC convergence diagnostic (Lambert & Vehtari 2022)
via
[`posterior::rstar()`](https://mc-stan.org/posterior/reference/rstar.html).
R\* measures whether a classifier can identify the chain that generated
a draw better than chance; values near 1 indicate convergence. Requires
per-chain posterior draws, which are extracted from the underlying fit
object (`fit_result$fit`):

- **BrmsFitter**: `fit_result$fit` is a `brmsfit`;
  [`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
  on it carries `.chain`.

- **CmdStanFitter**: `fit_result$fit$fit$draws()` returns a `draws_df`
  with `.chain`.

Fitters whose `fit_result$fit` carries no chain info (e.g.
`LinearRegressionFitter`, which has only a flat S x P draws matrix)
return `NA_real_` with a warning, since R\* is undefined without
multiple chains.

[`posterior::rstar()`](https://mc-stan.org/posterior/reference/rstar.html)
additionally requires the `caret` package (and a backend such as
`ranger` for the default random-forest classifier); on any error (e.g.
missing dependencies, too few chains) the metric degrades gracefully to
`NA_real_`.

Constructor for RstarMetric.

## Usage

``` r
RstarMetric(
  name = character(0),
  needs = character(0),
  required = FALSE,
  summary_type = "mean",
  uncertainty = FALSE,
  method = "rf"
)

rstar_metric(name = "rstar", uncertainty = FALSE, method = "rf")
```

## Arguments

- name:

  Character string naming the metric. Defaults to "rstar".

- uncertainty:

  Logical; see
  [`posterior::rstar()`](https://mc-stan.org/posterior/reference/rstar.html).
  Defaults to `FALSE`.

- method:

  Character; classifier passed to
  [`posterior::rstar()`](https://mc-stan.org/posterior/reference/rstar.html).

## Value

A `RstarMetric` object.

A `RstarMetric` object.

## Examples

``` r
rstar_metric()
#> <bayesim::RstarMetric>
#>  @ name        : chr "rstar"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ uncertainty : logi FALSE
#>  @ method      : chr "rf"
rstar_metric(uncertainty = FALSE)
#> <bayesim::RstarMetric>
#>  @ name        : chr "rstar"
#>  @ needs       : chr(0) 
#>  @ required    : logi FALSE
#>  @ summary_type: chr "mean"
#>  @ uncertainty : logi FALSE
#>  @ method      : chr "rf"
```
