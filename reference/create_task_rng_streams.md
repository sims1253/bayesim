# Create Task RNG Streams

Precomputes deterministic RNG streams for each task using L'Ecuyer-CMRG.
Each task receives its own independent RNG stream, ensuring
reproducibility regardless of execution order (sequential, parallel, or
resumed).

## Usage

``` r
create_task_rng_streams(global_seed, n_tasks)
```

## Arguments

- global_seed:

  Integer. The base seed for the simulation.

- n_tasks:

  Integer. Number of tasks to generate streams for.

## Value

A list of length `n_tasks`, where each element is an integer vector
representing the `.Random.seed` state for that task.

## Details

The caller's RNG kind and seed are restored on exit, including when
`.Random.seed` did not exist before the call.
