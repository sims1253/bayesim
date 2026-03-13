# Quantifies deviation from uniformity by the likelihood of observing the most extreme point on the empirical CDF of the given rank distribution according to Modrak et al. (2023, equation 7).

Modrak, Martin, Angie H. Moon, Shinyoung Kim, Paul Bürkner, Niko Huurre,
Kateřina Faltejsková, Andrew Gelman, and Aki Vehtari. "Simulation-Based
Calibration Checking for Bayesian Computation: The Choice of Test
Quantities Shapes Sensitivity." arXiv, June 15, 2023.
https://doi.org/10.48550/arXiv.2211.02383.

## Usage

``` r
gamma_discrepancy(ranks, post_warmup_draws, log = FALSE)
```

## Arguments

- ranks:

  Rank distribution

- post_warmup_draws:

  Number of posterior draws that were used to calculate the rank
  distribution.

- log:

  TRUE if the result should be on the log scale.

## Value

Measure quantifying deviation from uniformity. This value can be
compared to the distribution of gamma expected under uniformity
calculated by validation.gamma_null_distribution.

## Examples

``` r
# Compute gamma discrepancy for a sample of ranks
ranks <- sample(1:1000, 100, replace = TRUE)
gamma_discrepancy(ranks, post_warmup_draws = 1000)
#> [1] 0.02098557
```
