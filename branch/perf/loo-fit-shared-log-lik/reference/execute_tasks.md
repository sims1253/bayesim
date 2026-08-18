# Execute all tasks in grid

Execute all tasks in grid

## Usage

``` r
execute_tasks(
  task_grid,
  config,
  config_spec,
  fitter,
  metrics,
  retain,
  max_errors,
  progress,
  result_path = NULL,
  config_fingerprint = NULL,
  checkpoint_every = 50L,
  keep_checkpoints = 2L,
  prior_results_df = data.frame(task_id = character(), status = character()),
  prior_task_results = NULL,
  adaptive_next_check = NULL,
  adaptive_state = NULL,
  stop_on = NULL,
  verbose = TRUE,
  run_store = NULL
)
```

## Arguments

- task_grid:

  A task grid tibble from create_task_grid()

- config:

  A SimulationConfig S7 object

- config_spec:

  Plain list config spec for worker transport

- fitter:

  S7 Fitter object

- metrics:

  List of Metric objects

- retain:

  Character vector of what to retain

- max_errors:

  Maximum number of errors before stopping

- progress:

  Logical; if TRUE, show progress bar

- result_path:

  Character; path for checkpoint storage (optional)

- config_fingerprint:

  Character; configuration fingerprint for validation

- checkpoint_every:

  Integer; write checkpoint every N completed tasks. B4: also bounds the
  number of task results held in memory at once.

- keep_checkpoints:

  Integer; number of checkpoint commit directories retained. Pruning
  removes old commit directories only; immutable outcome shards and
  ledger history are never pruned.

- prior_results_df:

  Previously resumed result rows, cached in memory to avoid re-reading
  and re-hashing the prior checkpoint for every batch.

- stop_on:

  Optional adaptive stopping policy from the RunPolicy.

- verbose:

  Logical; if TRUE, print lifecycle messages independently of the task
  progress bar.

- run_store:

  Optional internal RunStore adapter used for checkpoint persistence.
  Defaults to an adapter created from `result_path`.

## Value

A list with task_results and task_grid
