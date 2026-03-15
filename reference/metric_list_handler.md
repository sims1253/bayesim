# Collect multiple metrics at once

A convenience function that collects all given metrics at once. This can
save time compared to manually calling all metric functions
individually, as some variables can be reused instead of being
calculated multiple times.

## Usage

``` r
metric_list_handler(fit, metrics, data_gen_output, fit_conf, ...)
```

## Arguments

- fit:

  A brmsfit object.

- metrics:

  A vector of metric identifiers. See details for supported identifiers.

- data_gen_output:

  Output from data generation containing true parameters

- fit_conf:

  Model configuration list or data.frame row

- ...:

  Additional arguments passed to metric calculation functions

## Value

A named list containing all requested metrics' results.

## Details

Currently, the following identifiers are supported:

- "bias": posterior bias

- "divergents": number of divergent transitions

- "ess": effective sample size

- "elpd_loo": ELPD-LOO

- "elpd_newdata": ELPD on new data

- "epred": expected predictions

- "mae_s": mean absolute error scaled

- "p_mean": posterior mean

- "p_sd": posterior standard deviation

- "pareto_k": Pareto k diagnostic values

- "pos_prob": posterior probability of positive values

- "ppred": posterior predictions

- "pq": posterior quantiles

- "q_true": true quantiles

- "r2_loo": LOO R-squared

- "r2_newdata": R-squared on new data

- "residuals": model residuals

- "rhat": R-hat diagnostic

- "rmse_loo": LOO root mean square error

- "rmse_newdata": RMSE on new data

- "rmse_s": scaled RMSE

- "rstar": R\* diagnostic

- "time_sampling": sampling time

- "time_total": total time

- "time_warmup": warmup time

- "y": observed values

Note,that not all identifiers are supported for each input class.

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires brms package and a fitted model
fit <- brms::brm(y ~ x, data = mydata)
metric_list_handler(fit, c("elpd_loo", "rhat"), data_gen_output, fit_conf)
} # }
```
