# Adjust gamma for simultaneous confidence intervals

Adjusts the coverage parameter to find simultaneous confidence intervals
for the ECDF of samples from the uniform distribution, as described in
Modrak et al. (2023).

## Usage

``` r
adjust_gamma(N, L, K = N, conf_level = 0.95)
```

## Source

This function is adapted from the SBC package
(https://github.com/hyunjimoon/SBC).

## Arguments

- N:

  Length of samples (chains).

- L:

  Number of samples (chains).

- K:

  Number of equally spaced evaluation points, i.e. the right ends of the
  partition intervals. Defaults to N.

- conf_level:

  Confidence level for the intervals. Default is 0.95.

## Value

The adjusted gamma value for computing confidence bands.

## Details

This function is used in Simulation-Based Calibration (SBC) to compute
confidence bands for rank statistics. It supports both single-chain
(L=1) and multi-chain analyses.

## References

Modrak, Martin, Angie H. Moon, Shinyoung Kim, Paul Bürkner, Niko Huurre,
Kateřina Faltejsková, Andrew Gelman, and Aki Vehtari. "Simulation-Based
Calibration Checking for Bayesian Computation: The Choice of Test
Quantities Shapes Sensitivity." arXiv, June 15, 2023.
https://doi.org/10.48550/arXiv.2211.02383.

## Examples

``` r
if (FALSE) { # \dontrun{
# Single chain
gamma <- adjust_gamma(N = 1000, L = 1)

# Multiple chains
gamma <- adjust_gamma(N = 250, L = 4)
} # }
```
