# Advance RNG stream

Advances the RNG state by n steps without returning random values.
Useful for skipping ahead in a stream.

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

This is a **pure function** with no side effects:

- It does NOT modify `.Random.seed` in the global environment

- It creates a local copy of the RNG state, advances it, and returns the
  new state

- The caller is responsible for setting the returned state if needed

## Examples

``` r
if (FALSE) { # \dontrun{
streams <- create_task_rng_streams(42, 10)
advanced <- advance_rng_stream(streams[[1]], n = 5)
} # }
```
