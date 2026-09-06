# Build final simulation result

Build final simulation result

## Usage

``` r
build_simulation_result(
  config,
  task_results,
  task_grid,
  final_results_df = NULL,
  elapsed,
  checkpoint_path = NULL
)
```

## Arguments

- config:

  A SimulationConfig S7 object

- task_results:

  List of bayesim_task_result objects

- task_grid:

  The task grid tibble with updated statuses

- final_results_df:

  Dataframe with merged results (prior + new)

- elapsed:

  Total elapsed time in seconds

- checkpoint_path:

  Character; path where checkpoints were stored

## Value

A bayesim_simulation_result S3 object
