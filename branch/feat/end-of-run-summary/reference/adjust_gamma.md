# Adjust the coverage parameter for simultaneous ECDF confidence bands

Computes the gamma coverage parameter such that the simultaneous
confidence envelope of the ECDF of a uniform sample of size N has
(approximately) the requested confidence level (Säilynoja et al. 2022).
For L = 1 (single chain/sample) the result is exact via dynamic
programming; for L \> 1 the 0.x code used a Monte-Carlo simulation that
depended on an external `u_scale` helper not ported here, so the exact L
= 1 band is used (slightly conservative).

## Usage

``` r
adjust_gamma(N, L, K = N, conf_level = 0.95)
```

## Arguments

- N:

  Integer; number of samples (ranks).

- L:

  Integer; number of samples/chains. Default 1.

- K:

  Integer; number of equally spaced evaluation points (right ends of the
  partition intervals). Defaults to N.

- conf_level:

  Numeric in (0,1); confidence level. Default 0.95.

## Value

Numeric gamma in (0, 1 - conf_level).
