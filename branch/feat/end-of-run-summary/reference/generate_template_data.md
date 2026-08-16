# Generate template data for prefit compilation

Produces ONE small dataset using the configured data generator with a
throwaway RNG stream. The dataset is discarded after compilation; its
only purpose is to fix the Stan data structure (variable names and
types) so the compiled binary matches the real task data.

## Usage

``` r
generate_template_data(data_generator, data_spec)
```

## Arguments

- data_generator:

  The simulation data generator function.

- data_spec:

  A named list (one row of `data_grid` as a list).

## Value

A data.frame suitable for `brms::brm(data =)`.

## Details

Calls `data_generator(data_spec, task_ctx)` and returns its `$train`
data frame. `task_ctx$template` is TRUE so generators can recognize this
structural template call.
