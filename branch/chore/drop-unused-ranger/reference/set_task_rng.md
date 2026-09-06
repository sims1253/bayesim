# Set RNG state for a task

Restores .Random.seed from a precomputed stream. Call this at the start
of each worker task.

## Usage

``` r
set_task_rng(rng_stream)
```

## Arguments

- rng_stream:

  Integer vector from create_task_rng_streams()

## Value

NULL (invisibly). Side effect: sets `.Random.seed` in global
environment.

## Examples

``` r
if (FALSE) { # \dontrun{
streams <- create_task_rng_streams(42, 10)
set_task_rng(streams[[1]])
} # }
```
