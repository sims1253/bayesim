# Execute a single simulation task

Runs one task: generate data, fit model, compute metrics. All errors are
captured and converted to task results.

## Usage

``` r
run_task(
  task,
  config_spec,
  fitter,
  metrics,
  retain = c("metrics", "diagnostics")
)
```

## Arguments

- task:

  A task specification list (from get_task_spec) containing:

  - `task_id`: Unique task identifier string

  - `data_spec`: Data generation specification

  - `fit_spec`: Model fitting specification

  - `rng_seed`: Precomputed RNG seed for this task

  - `task_ctx`: Task context with indices (data_idx, fit_idx, rep_idx)

- config_spec:

  Plain list config spec (for worker transport) containing:

  - `data_generator`: Function to generate data

  - Other configuration parameters

- fitter:

  S7 Fitter object that implements the fit_model(), predict_fit(),
  log_lik_matrix(), and loo_fit() methods

- metrics:

  List of S7 Metric objects to compute

- retain:

  Character vector of what to retain in results. Options:

  - "metrics": Always retained (default)

  - "diagnostics": Fit diagnostics (default)

  - "fit": The raw fit object

  - "draws": Posterior draws matrix

## Value

A bayesim_task_result S3 object containing:

- `task_id`: Task identifier

- `status`: "success" or "failed"

- `metrics`: Named list of computed metrics (if successful)

- `diagnostics`: Named list of fit diagnostics (if successful)

- `timing`: List with `total` elapsed time

- `warnings`: Character vector of warning messages

- `error`: Error information (if failed)

## Details

The task execution follows these steps:

1.  Set the RNG state for deterministic reproducibility

2.  Generate data using the data_generator function

3.  Validate the data bundle

4.  Fit the model using the fitter

5.  Build context for metric computation (predictions, log_lik, loo)

6.  Compute all metrics

7.  Apply retention policy to manage memory

Error handling:

- Data generation errors are normalized to bayesim_data_error

- Fit errors are normalized to bayesim_fit_error

- Metric errors are handled in compute_all_metrics()

- Fatal errors (config, contract) propagate and stop the simulation

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a task specification
task <- list(
  task_id = "d001_f001_r00001",
  data_spec = list(n = 100, effect_size = 0.5),
  fit_spec = list(prior = "normal(0, 1)"),
  rng_seed = create_task_rng_streams(42, 1)[[1]],
  task_ctx = list(data_idx = 1, fit_idx = 1, rep_idx = 1)
)

# Run the task
result <- run_task(
  task = task,
  config_spec = config_spec,
  fitter = my_fitter,
  metrics = list(rmse_metric, coverage_metric),
  retain = c("metrics", "diagnostics")
)

# Check result
if (result$status == "success") {
  print(result$metrics)
} else {
print(result$error$error_message)
}
} # }
```
