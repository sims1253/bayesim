# Required number of replicates for a target MCSE

Inverts the Monte-Carlo standard error (MCSE) formulas used by
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)
to return the number of replicates `n` required to achieve a target
MCSE. Useful for planning a simulation study (a precision / power
calculation) before running it.

The MCSE formulas (Morris et al. 2019; rsimsum) invert to:

- **coverage** (binary 0/1 outcome): `MCSE = sqrt(p (1 - p) / n)` =\>
  `n = p (1 - p) / MCSE^2`. The default `p = 0.5` is the conservative
  max-variance case (it maximises `p (1 - p)` and therefore gives the
  largest required `n`).

- **continuous** metrics (bias / mean / model SE): `MCSE = sd / sqrt(n)`
  =\> `n = (sd / MCSE)^2`. Requires an assumed standard deviation
  `assumed_sd` for the per-replicate point estimate (e.g. a guess at the
  empirical SE of the estimator).

The returned value is `ceiling(n)` so it is always a whole number of
replicates.

## Usage

``` r
n_replicates_for_target(
  target_mcse,
  metric_type = c("coverage", "continuous"),
  p = 0.5,
  assumed_sd = NULL
)
```

## Arguments

- target_mcse:

  Numeric scalar \> 0. The MCSE you want to achieve.

- metric_type:

  Character scalar: `"coverage"` (binary coverage metrics) or
  `"continuous"` (bias / mean / model SE metrics).

- p:

  Numeric scalar in `[0, 1]`. Assumed coverage probability. Only used
  when `metric_type = "coverage"`. Defaults to `0.5`, the
  variance-maximising (most conservative) choice.

- assumed_sd:

  Numeric scalar \> 0. Assumed standard deviation of the per-replicate
  point estimate. Required when `metric_type = "continuous"`; ignored
  otherwise.

## Value

An integer scalar: the number of replicates required (the ceiling of the
inverted-MCSE `n`).

## See also

[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)
for the MCSE formulas being inverted.

## Examples

``` r
if (FALSE) { # \dontrun{
# Coverage: target MCSE 0.03 under the conservative p = 0.5.
n_replicates_for_target(0.03, "coverage")

# Coverage at an assumed rate of 0.9.
n_replicates_for_target(0.03, "coverage", p = 0.9)

# Continuous metric (e.g. bias) with assumed sd of the point estimate.
n_replicates_for_target(0.05, "continuous", assumed_sd = 0.5)
} # }
```
