# Validate the optional adaptive-stopping policy (I3)

`NULL` is valid (no adaptive stopping). Otherwise the input must be a
list with: `estimand` (character), `measure` (one of
`VALID_STOP_MEASURES`), `target_mcse` (numeric \> 0), and optional
`min_reps` (integer, default 50) and `check_every` (integer, default
50). Returns a normalized list.

## Usage

``` r
validate_stop_on(stop_on)
```

## Arguments

- stop_on:

  NULL or a list.

## Value

NULL or a normalized list.
