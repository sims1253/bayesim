# Create a Simulation Configuration

Creates a SimulationConfig S7 object that fully specifies a simulation
study. This configuration defines the data generation grid, fitting
grid, metrics, and execution parameters.

## Usage

``` r
simulation_config(
  data_grid = NULL,
  fit_grid = NULL,
  task_grid = NULL,
  data_generator,
  fitter = NULL,
  metrics = NULL,
  n_replicates = 1L,
  seed,
  result_path = NULL,
  checkpoint_format = c("rds", "parquet"),
  checkpoint_every = 50L,
  chunk_size = NULL,
  max_in_memory = NULL,
  retain = c("metrics", "diagnostics"),
  max_errors = Inf
)
```

## Arguments

- data_grid:

  A data.frame with data generation specifications. Each row represents
  a distinct data configuration to simulate.

- fit_grid:

  A data.frame with model fitting specifications. Each row represents a
  distinct model configuration to fit.

- task_grid:

  Optional pre-computed task grid. If provided, overrides data_grid and
  fit_grid. Must contain either data_spec/fit_spec list-columns or
  data_idx/fit_idx index columns.

- data_generator:

  A function with signature
  `(data_spec, seed, task_ctx) -> data_bundle`. Generates data for a
  single replicate given a data specification row.

- fitter:

  An S7 Fitter object that handles model fitting.

- metrics:

  A list of Metric objects.

- n_replicates:

  Positive integer. Number of replicates per data/fit combination.

- seed:

  Integer. Base seed for reproducible random number generation.

- result_path:

  NULL or character path. If provided, results are saved here.

- checkpoint_format:

  Character scalar. Checkpoint storage format. Currently only `"rds"` is
  implemented for checkpoint persistence.

- checkpoint_every:

  Positive integer. Save progress every N tasks.

- chunk_size:

  Positive integer. Maximum number of task results to keep in memory
  before forcing a checkpoint write. Defaults to `checkpoint_every`.

- max_in_memory:

  Deprecated alias for `chunk_size`.

- retain:

  Character vector. What to retain in results. Must be subset of
  `c("metrics", "diagnostics", "draws", "predictions", "fit", "data", "warnings")`.

- max_errors:

  Numeric. Maximum errors before stopping. Use `Inf` for no limit.

## Value

An S7 SimulationConfig object.

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(
  data_grid = data.frame(n = c(100, 500), effect = c(0.5, 1.0)),
  fit_grid = data.frame(model = c("baseline", "full")),
  data_generator = my_data_gen,
  fitter = my_fitter,
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 100L,
  seed = 42L,
  checkpoint_format = "rds"
)
} # }
```
