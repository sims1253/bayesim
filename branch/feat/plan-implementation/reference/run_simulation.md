# Run a Simulation Study

Executes a complete simulation study with deterministic reproducibility.

## Usage

``` r
run_simulation(
  config,
  resume = c("auto", "never", "must"),
  progress = TRUE,
  workers = NULL,
  verbose = TRUE
)
```

## Arguments

- config:

  A SimulationConfig S7 object

- resume:

  Character strategy controlling how an existing `result_path` is
  treated: `"auto"` (default) resumes when the path holds a compatible
  run and starts fresh only when the path is absent or empty — an
  existing run that is incompatible or corrupt aborts rather than being
  silently overwritten, as does a non-empty directory without a run
  manifest; `"never"` starts fresh and errors if `result_path` already
  holds a run or unrelated files; `"must"` resumes and errors when no
  compatible checkpoint exists. Only tasks with terminal status
  (`"success"`/`"failed"`) are carried over; all other tasks re-run with
  their original RNG streams.

- progress:

  Logical; if TRUE, show progress bar

- workers:

  Positive integer, NULL, or "multisession". When non-NULL,
  `mirai::daemons(workers)` is set up for the run and torn down on exit
  — the simple path for local parallelism. `workers = 1` is genuinely
  sequential: no daemons are launched, so package-external S7 fitters
  and metrics keep their method dispatch (S7 method tables are
  registered per process and do not travel to daemon workers). Parallel
  execution starts at `workers >= 2`. Must be NULL when daemons are
  already set (use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  directly for the advanced/HPC path: remote daemons, TLS, etc.).
  Daemons are set before the model bank ships.

- verbose:

  Logical; if TRUE (default), print preflight, lifecycle, and completion
  messages. This is independent of `progress`: set `progress = FALSE` to
  hide the task bar while keeping run messages, or `verbose = FALSE` for
  a quiet programmatic run.

## Value

A bayesim_simulation_result S3 object

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(
  data_grid = data.frame(n = c(100, 500)),
  fit_grid = data.frame(model = "baseline"),
  data_generator = my_data_gen,
  fitter = my_fitter,
  metrics = list(pred_rmse_metric(), pred_bias_metric()),
  n_replicates = 100L,
  seed = 42L
)
result <- run_simulation(config)
} # }
```
