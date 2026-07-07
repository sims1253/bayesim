# Generate template data for prefit compilation

Produces ONE small dataset using the configured data generator with a
throwaway RNG stream. The dataset is discarded after compilation; its
only purpose is to fix the Stan data structure (variable names and
types) so the compiled binary matches the real task data.

## Usage

``` r
generate_template_data(data_generator, data_spec, seed = 0L)
```

## Arguments

- data_generator:

  The simulation data generator function.

- data_spec:

  A named list (one row of `data_grid` as a list).

- seed:

  Integer throwaway seed (default 0L).

## Value

A data.frame suitable for `brms::brm(data =)`.

## Details

Calls `data_generator(data_spec, seed, task_ctx)` and returns its
`$train` data frame. The seed is a fixed throwaway constant (NOT the
simulation seed) so the same template data is produced on controller and
daemons.
