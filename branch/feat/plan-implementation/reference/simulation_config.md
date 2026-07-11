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
  checkpoint_format = c("rds"),
  checkpoint_every = 50L,
  keep_checkpoints = 2L,
  retain = c("metrics", "diagnostics"),
  max_errors = Inf,
  daemon_setup = NULL,
  stop_on = NULL,
  summary_format = c("rds", "parquet")
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

  A function with signature `(data_spec, task_ctx) -> data_bundle`.
  Generates data for a single replicate given a data specification row.
  `task_ctx$seed` carries the per-task integer seed for backends that
  need one.

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
  implemented for checkpoint persistence. (B4: excluded from the config
  fingerprint — it is runtime policy.)

- checkpoint_every:

  Positive integer. Save progress every N tasks. This single knob also
  bounds the number of task results held in memory at once (B4: the
  former separate `chunk_size` knob was merged into this).

- keep_checkpoints:

  Positive integer. Number of complete checkpoint snapshots to retain on
  disk. Defaults to 2, preserving the latest snapshot plus one fallback
  for corruption recovery. Runtime policy; excluded from the config
  fingerprint.

- retain:

  Character vector. What to retain in results. Must be subset of
  `c("metrics", "diagnostics", "draws", "predictions", "fit", "data", "warnings")`.
  (B4: excluded from the config fingerprint — changing retention must
  not invalidate resume.)

- max_errors:

  Numeric. Maximum errors before stopping. Use `Inf` for no limit. (B4:
  excluded from the config fingerprint.)

- daemon_setup:

  Optional function run once per mirai daemon (via
  [`mirai::everywhere()`](https://mirai.r-lib.org/reference/everywhere.html))
  before tasks start, e.g. to configure cmdstan paths or load a model
  bank. Ignored when no daemons are set. Default NULL.

- stop_on:

  Optional adaptive stopping policy (experimental). `NULL` (default)
  runs all tasks. Otherwise a list with elements: `estimand` (character
  parameter name), `measure` (one of `"bias"`, `"coverage"`, `"emp_se"`,
  `"mse"`, `"model_se"`), `target_mcse` (numeric \> 0), `min_reps`
  (integer, default 50), `check_every` (integer, default 50). Once the
  MCSE of `measure` for `estimand` falls below `target_mcse` AND at
  least `min_reps` replicates have completed, remaining pending tasks
  are marked `"skipped"` and the run stops. (I3: excluded from the
  config fingerprint — it is runtime policy.)

- summary_format:

  Character scalar. Output format for the final summary. `"rds"`
  (default) writes nothing extra – the canonical rds checkpoint carries
  the summary and remains the resume artifact. `"parquet"` additionally
  writes `<result_path>/summary.parquet` using the suggested
  `nanoparquet` package, for downstream consumption (pandas, arrow,
  polars). (I8: excluded from the config fingerprint – runtime policy.)

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
  metrics = list(pred_rmse_metric(), pred_bias_metric()),
  n_replicates = 100L,
  seed = 42L,
  checkpoint_format = "rds"
)
} # }
```
