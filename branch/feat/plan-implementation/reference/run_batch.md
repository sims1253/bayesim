# Run one batch of tasks

Dispatches a batch of simulation tasks via purrr's mirai integration
([`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) +
[`purrr::in_parallel()`](https://purrr.tidyverse.org/reference/in_parallel.html)).
With no daemons set, purrr falls back to sequential execution
automatically, so there is a single code path (C1). mirai remains the
daemon engine; daemons/model bank/`daemon_setup` are managed by
[`execute_tasks()`](https://sims1253.github.io/bayesim/reference/execute_tasks.md).

## Usage

``` r
run_batch(batch_tasks, config_spec, fitter, metrics, retain, progress = FALSE)
```

## Arguments

- batch_tasks:

  List of task specification lists

- config_spec:

  Plain list config spec for worker transport

- fitter:

  S7 Fitter object

- metrics:

  List of Metric objects

- retain:

  Character vector of what to retain

- progress:

  Logical; if TRUE, surface purrr's progress display

## Value

A list of bayesim_task_result objects

## Details

[`run_task_safe()`](https://sims1253.github.io/bayesim/reference/run_task_safe.md)
is total (never throws), so transport is pure transport: every returned
element is a `bayesim_task_result`. Transport-level failures (e.g.
daemon death) surface as errors from
[`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) and are
re-raised as a bayesim fatal error by the caller.
