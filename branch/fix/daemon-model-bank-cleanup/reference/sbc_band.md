# Simultaneous confidence band for a uniform ECDF.

Returns the lower/upper bounds of the simultaneous ECDF confidence
envelope (Säilynoja et al. 2022) at K equally spaced evaluation points,
given N samples and a confidence level.

## Usage

``` r
sbc_band(N, K = N, conf_level = 0.95)
```

## Arguments

- N:

  Integer; number of samples (ranks).

- K:

  Integer; number of evaluation points. Defaults to N.

- conf_level:

  Numeric in (0,1); confidence level.

## Value

A list with `x` (the grid 0:K / K) and `lower` and `upper` numeric
vectors of length K + 1 over that grid.
