# Lighten task result after checkpointing

Removes heavy objects from a task result that has been checkpointed.
This reduces memory usage while keeping lightweight summary data needed
for the final result construction.

## Usage

``` r
lighten_task_result(task_result, retain)
```

## Arguments

- task_result:

  A bayesim_task_result object

- retain:

  Character vector of retention options specifying what to keep

## Value

A lightened bayesim_task_result object with heavy fields removed. If all
heavy fields are removed and only minimal data remains, may return a
lightweight summary object instead.

## Details

After a task result has been checkpointed to disk, we can safely remove
heavy objects from memory since they can be reconstructed from the
checkpoint. This function:

- Always keeps: task_id, status, metrics, timing, error

- Keeps diagnostics only if "diagnostics" is in retain

- Keeps retained fit, draws, predictions, and data fields; these are
  intentionally omitted only when their retention option is absent

- Keeps warnings only if "warnings" is in retain

The retention policy is authoritative: metrics-only runs stay
lightweight, while users who explicitly retain heavy artifacts receive
them in the final result even when checkpointing is enabled.

## See also

[`execute_tasks()`](https://sims1253.github.io/bayesim/reference/execute_tasks.md),
[`write_checkpoint()`](https://sims1253.github.io/bayesim/reference/write_checkpoint.md)
