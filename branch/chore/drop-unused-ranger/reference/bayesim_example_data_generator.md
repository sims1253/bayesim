# Built-in example data generator

Generates a simple linear-regression dataset (`y = beta * x + noise`).

## Usage

``` r
bayesim_example_data_generator(data_spec, task_ctx)
```

## Arguments

- data_spec:

  List with `n`, `beta`, `sigma`.

- task_ctx:

  Task context (carries `seed` for backends that need one).

## Value

A `data_bundle` list.

## Details

Generators consume the ambient RNG state (the worker restores the
per-task L'Ecuyer stream via
[`set_task_rng()`](https://sims1253.github.io/bayesim/reference/set_task_rng.md)).
`task_ctx$seed` carries an integer seed for backends that need one;
generators must not use it to re-seed.
