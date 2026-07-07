# Advance RNG stream

Advances the RNG state by n steps and returns the new state.

## Usage

``` r
advance_rng_stream(rng_stream, n = 1L)
```

## Arguments

- rng_stream:

  Integer vector RNG state

- n:

  Number of steps to advance (default: 1)

## Value

Advanced RNG state (integer vector)

## Details

This function temporarily sets `.Random.seed` in `.GlobalEnv` to advance
the stream, then returns the new state. It does NOT permanently modify
the global `.Random.seed` (unless it did not exist before, in which case
it is removed after use). The caller is responsible for applying the
returned state if needed (e.g., via
[`set_task_rng()`](https://sims1253.github.io/bayesim/reference/set_task_rng.md)).

## Examples

``` r
if (FALSE) { # \dontrun{
streams <- create_task_rng_streams(42, 10)
advanced <- advance_rng_stream(streams[[1]], n = 5)
} # }
```
