# Full forward sampling of a brms model's responses for a SINGLE posterior draw

Simulates the response(s) of a brmsfit at a single posterior draw,
respecting dependency order among multiple responses. Ported from the
SBC package and the 0.x bayesim codebase; bayeshear/SBC/future
dependencies removed.

## Usage

``` r
brms_full_ppred(fit, newdata = NULL, draw = 1L)
```

## Arguments

- fit:

  A brmsfit with posterior draws.

- newdata:

  Optional data.frame of predictors. If `NULL`, uses `fit$data`.

- draw:

  Single integer draw index to simulate.

## Value

A data.frame (a copy of `newdata`) with the simulated response column(s)
filled in.

## Details

[`ifs_generator()`](https://sims1253.github.io/bayesim/reference/ifs_generator.md)
only ever needs one draw per task (the deterministic `rep_idx`-indexed
draw), so this function takes a single `draw` index and returns a single
data.frame. This eliminates the prior index mismatch where results were
stored at `pp_data[[draw_value]]` but read at `simulated[[1]]`.
