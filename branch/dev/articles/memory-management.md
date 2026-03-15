# Memory Management for Large Simulations

## Introduction

Simulation studies in Bayesian modeling can consume large amounts of
memory. When running thousands of tasks, each producing posterior draws,
fitted model objects, predictions, and diagnostics, the accumulated data
can exhaust available RAM. This vignette explains how to use bayesim’s
retention system to control memory usage during large-scale simulations.

Consider a simulation with: - 50 data configurations - 5 model
variants - 100 replicates per combination - 4,000 posterior draws per
fit - 50 parameters tracked

This creates 25,000 tasks. If each fit object is 50 MB and each draws
matrix is 20 MB, retaining everything would require over 1.7 TB of
memory. We need a way to retain only what is necessary.

## The Retention System

### What Gets Stored

During simulation, each task produces several types of artifacts:

| Artifact      | Description                                         | Typical Size  |
|---------------|-----------------------------------------------------|---------------|
| `metrics`     | Computed metric values (RMSE, bias, coverage, etc.) | Small (KB)    |
| `diagnostics` | Convergence diagnostics (R-hat, ESS, divergences)   | Small (KB)    |
| `draws`       | Posterior draws matrix                              | Medium (MB)   |
| `predictions` | Predicted values for test/training data             | Variable      |
| `fit`         | Raw fit object from the modeling backend            | Large (MB-GB) |
| `data`        | Input training and test data                        | Variable      |
| `warnings`    | Warning messages from fitting                       | Small (KB)    |

### The `retain` Parameter

The `retain` parameter in
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
controls which artifacts are kept in memory:

``` r
library(bayesim)

# Default retention: metrics and diagnostics
config <- simulation_config(
  data_grid = data.frame(n = 100, sigma = 1),
  fit_grid = data.frame(model = "linear"),
  data_generator = function(data_spec, seed, task_ctx) {
    # Note: seed is a scalar task seed and the task RNG stream is restored
    # before each call.
    n <- data_spec$n
    x <- rnorm(n)
    y <- 1 + 2 * x + rnorm(n, sd = data_spec$sigma)
    list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(intercept = 1, slope = 2, sigma = data_spec$sigma),
      vars_of_interest = c("intercept", "slope", "sigma"),
      references = c(intercept = 0, slope = 0, sigma = 1),
      meta = list()
    )
  },
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 10L,
  seed = 42L,
  retain = c("metrics", "diagnostics")  # Default
)
```

### Retention Profiles

bayesim provides three common retention profiles for typical use cases.
You can use the
[`resolve_retention()`](https://sims1253.github.io/bayesim/reference/resolve_retention.md)
function to get the explicit character vector for each profile:

#### Minimal Profile: Keep Only Essential Metrics

Use this for very large simulations where you only need summary
statistics:

``` r
# Get the minimal retention profile
minimal_retain <- resolve_retention("minimal")
print(minimal_retain)
#> [1] "metrics"

config_minimal <- simulation_config(
  data_grid = data.frame(n = 100, sigma = 1),
  fit_grid = data.frame(model = "linear"),
  data_generator = function(data_spec, seed, task_ctx) {
    # Note: seed is a scalar task seed and the task RNG stream is restored
    # before each call.
    n <- data_spec$n
    x <- rnorm(n)
    y <- 1 + 2 * x + rnorm(n, sd = data_spec$sigma)
    list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(intercept = 1, slope = 2, sigma = data_spec$sigma),
      vars_of_interest = c("intercept", "slope", "sigma"),
      references = c(intercept = 0, slope = 0, sigma = 1),
      meta = list()
    )
  },
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 10L,
  seed = 42L,
  retain = minimal_retain  # Only metrics
)
```

The `"minimal"` profile retains only `metrics`, giving you the smallest
memory footprint.

#### Standard Profile: Keep Metrics, Diagnostics, and Warnings

This is the default profile, balancing information with memory usage:

``` r
# Get the standard retention profile
standard_retain <- resolve_retention("standard")
print(standard_retain)
#> [1] "metrics"     "diagnostics" "warnings"

config_standard <- simulation_config(
  data_grid = data.frame(n = 100, sigma = 1),
  fit_grid = data.frame(model = "linear"),
  data_generator = function(data_spec, seed, task_ctx) {
    # Note: seed is a scalar task seed and the task RNG stream is restored
    # before each call.
    n <- data_spec$n
    x <- rnorm(n)
    y <- 1 + 2 * x + rnorm(n, sd = data_spec$sigma)
    list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(intercept = 1, slope = 2, sigma = data_spec$sigma),
      vars_of_interest = c("intercept", "slope", "sigma"),
      references = c(intercept = 0, slope = 0, sigma = 1),
      meta = list()
    )
  },
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 10L,
  seed = 42L,
  retain = standard_retain
)
```

The `"standard"` profile retains `metrics`, `diagnostics`, and
`warnings`, useful for checking convergence and capturing any warning
messages across your simulation.

#### Debug Profile: Keep Everything

Use this for debugging or small studies where you need full access to
all artifacts:

``` r
# Get the debug retention profile
debug_retain <- resolve_retention("debug")
print(debug_retain)
#> [1] "metrics"     "diagnostics" "draws"       "predictions" "fit"        
#> [6] "data"        "warnings"

config_debug <- simulation_config(
  data_grid = data.frame(n = 100, sigma = 1),
  fit_grid = data.frame(model = "linear"),
  data_generator = function(data_spec, seed, task_ctx) {
    # Note: seed is a scalar task seed and the task RNG stream is restored
    # before each call.
    n <- data_spec$n
    x <- rnorm(n)
    y <- 1 + 2 * x + rnorm(n, sd = data_spec$sigma)
    list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(intercept = 1, slope = 2, sigma = data_spec$sigma),
      vars_of_interest = c("intercept", "slope", "sigma"),
      references = c(intercept = 0, slope = 0, sigma = 1),
      meta = list()
    )
  },
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 10L,
  seed = 42L,
  retain = debug_retain  # Everything
)
```

The `"debug"` profile retains all artifacts: `metrics`, `diagnostics`,
`draws`, `predictions`, `fit`, `data`, and `warnings`.

## Custom Retention

For fine-grained control, specify a custom character vector:

``` r
# Keep metrics, diagnostics, and draws, but not fit objects or data
config_custom <- simulation_config(
  data_grid = data.frame(n = 100, sigma = 1),
  fit_grid = data.frame(model = "linear"),
  data_generator = function(data_spec, seed, task_ctx) {
    # Note: seed is a scalar task seed and the task RNG stream is restored
    # before each call.
    n <- data_spec$n
    x <- rnorm(n)
    y <- 1 + 2 * x + rnorm(n, sd = data_spec$sigma)
    list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(intercept = 1, slope = 2, sigma = data_spec$sigma),
      vars_of_interest = c("intercept", "slope", "sigma"),
      references = c(intercept = 0, slope = 0, sigma = 1),
      meta = list()
    )
  },
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric(), coverage_metric()),
  n_replicates = 10L,
  seed = 42L,
  retain = c("metrics", "diagnostics", "draws", "warnings")
)
```

Valid retention options are: - `"metrics"` - Always retained
(required) - `"diagnostics"` - Convergence diagnostics - `"draws"` -
Posterior draws matrix - `"predictions"` - Predicted values - `"fit"` -
Raw fit object - `"data"` - Input data - `"warnings"` - Warning messages

## When to Use What

### Small Studies (\< 100 Tasks)

For small studies, you can safely use the debug profile:

``` r
# Small study: 2 data configs × 1 model × 20 replicates = 40 tasks
small_config <- simulation_config(
  data_grid = data.frame(n = c(100, 500), sigma = 1),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_data_generator,
  fitter = my_fitter,
  metrics = list(rmse_metric(), coverage_metric()),
  n_replicates = 20L,
  seed = 42L,
  retain = resolve_retention("debug")  # Keep everything for detailed analysis
)
```

With fewer than 100 tasks, even retaining all artifacts typically uses
less than 5 GB of memory.

### Medium Studies (100-1000 Tasks)

For medium-sized studies, use the standard profile to keep diagnostics:

``` r
# Medium study: 10 data configs × 2 models × 50 replicates = 1,000 tasks
medium_config <- simulation_config(
  data_grid = data.frame(
    n = rep(c(100, 200, 500), each = 2),
    sigma = rep(c(0.5, 1), 3)
  ),
  fit_grid = data.frame(model = c("baseline", "full")),
  data_generator = my_data_generator,
  fitter = my_fitter,
  metrics = list(rmse_metric(), bias_metric(), coverage_metric()),
  n_replicates = 50L,
  seed = 42L,
  retain = resolve_retention("standard")
)
```

The `"standard"` profile keeps enough information to monitor convergence
while avoiding the memory overhead of storing all fit objects.

### Large Studies (\> 1000 Tasks)

For large studies, use the minimal profile or a custom profile:

``` r
# Large study: 20 data configs × 5 models × 100 replicates = 10,000 tasks
large_config <- simulation_config(
  data_grid = expand.grid(
    n = c(50, 100, 200, 500),
    sigma = c(0.5, 1, 2),
    effect = c(0.5, 1.0)
  ),
  fit_grid = data.frame(model = c("m1", "m2", "m3", "m4", "m5")),
  data_generator = my_data_generator,
  fitter = my_fitter,
  metrics = list(rmse_metric(), coverage_metric()),
  n_replicates = 100L,
  seed = 42L,
  retain = resolve_retention("minimal"),  # Only metrics
  result_path = "large_simulation",
  checkpoint_every = 100L
)
```

With `"minimal"`, you keep only the computed metrics, reducing memory
usage to megabytes even for thousands of tasks.

## Checkpointing for Large Studies

When running large simulations, combine retention profiles with
checkpointing to manage both memory and execution time:

``` r
# Very large study with checkpointing
huge_config <- simulation_config(
  data_grid = expand.grid(
    n = c(100, 500, 1000),
    sigma = c(0.5, 1, 2),
    effect = c(0.5, 1.0, 2.0)
  ),
  fit_grid = data.frame(model = c("baseline", "complex", "spline")),
  data_generator = my_data_generator,
  fitter = my_fitter,
  metrics = list(rmse_metric(), bias_metric(), coverage_metric()),
  n_replicates = 200L,
  seed = 42L,
  retain = resolve_retention("minimal"),
  result_path = "huge_simulation_results",  # Enable checkpointing
  checkpoint_every = 50L  # Save every 50 tasks
)

# Run with resume capability
result <- run_simulation(huge_config, resume = "auto")
```

### How Checkpointing Manages Memory

When checkpointing is enabled, bayesim implements memory-bounded
execution:

1.  Task results are accumulated in memory during execution
2.  When `chunk_size` tasks complete (or `checkpoint_every`, whichever
    comes first), a checkpoint is written to disk
3.  After checkpointing, heavy objects in task results are cleared from
    memory while keeping lightweight summary data (metrics, diagnostics)
4.  Only pending tasks and lightweight summaries remain in memory

Memory usage stays bounded by `chunk_size` tasks, regardless of total
study size.

You can control the memory bound with the `chunk_size` parameter:

``` r
# Tighter memory control - checkpoint and clear every 25 tasks
tight_memory_config <- simulation_config(
  data_grid = expand.grid(n = c(100, 500), sigma = c(0.5, 1)),
  fit_grid = data.frame(model = c("m1", "m2")),
  data_generator = my_data_generator,
  fitter = my_fitter,
  metrics = list(rmse_metric(), coverage_metric()),
  n_replicates = 100L,
  seed = 42L,
  retain = resolve_retention("minimal"),
  result_path = "tight_memory_simulation",
  checkpoint_every = 100L,  # Write checkpoint every 100 tasks
  chunk_size = 25L          # Clear memory after 25 tasks
)
```

When `chunk_size` is smaller than `checkpoint_every`, the system
checkpoints early (at the `chunk_size` threshold) to free memory. When
`chunk_size` is larger, checkpoints happen at `checkpoint_every`
intervals as usual. `max_in_memory` remains available as a compatibility
alias.

### Resuming Efficiently

If a simulation is interrupted, resuming continues from the last
checkpoint without reloading previous full results:

``` r
# First run - gets interrupted
result <- run_simulation(huge_config, resume = "auto")

# Later - resume from checkpoint
result <- resume_simulation("huge_simulation_results")
```

The checkpoint system stores summary data (metrics, diagnostics) needed
for analysis. Heavy high-cardinality metric payloads are externalized
into artifact files under the result directory rather than widening the
main summary table indefinitely.

## Example: Comparing Retention Profiles

Let’s compare memory usage across different retention profiles:

``` r
# Data generator for comparison
demo_generator <- function(data_spec, seed, task_ctx) {
  # Note: seed is a scalar task seed and the task RNG stream is restored
  # before each call.
  n <- data_spec$n
  x <- rnorm(n)
  y <- data_spec$intercept + data_spec$slope * x + rnorm(n, sd = data_spec$sigma)
  
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(
      intercept = data_spec$intercept,
      slope = data_spec$slope,
      sigma = data_spec$sigma
    ),
    vars_of_interest = c("intercept", "slope", "sigma"),
    references = c(intercept = 0, slope = 0, sigma = 1),
    meta = list()
  )
}

# Run with minimal retention
config_minimal <- simulation_config(
  data_grid = data.frame(
    n = 100,
    intercept = 1,
    slope = 2,
    sigma = 1
  ),
  fit_grid = data.frame(model = "linear"),
  data_generator = demo_generator,
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric(), coverage_metric()),
  n_replicates = 20L,
  seed = 42L,
  retain = resolve_retention("minimal")
)

result_minimal <- run_simulation(config_minimal, progress = FALSE)
#> ℹ Starting simulation with 20 tasks

# Run with standard retention
config_standard <- simulation_config(
  data_grid = data.frame(
    n = 100,
    intercept = 1,
    slope = 2,
    sigma = 1
  ),
  fit_grid = data.frame(model = "linear"),
  data_generator = demo_generator,
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric(), coverage_metric()),
  n_replicates = 20L,
  seed = 42L,
  retain = resolve_retention("standard")
)

result_standard <- run_simulation(config_standard, progress = FALSE)
#> ℹ Starting simulation with 20 tasks

# Compare summary outputs
cat("=== Minimal Retention Summary ===\n")
#> === Minimal Retention Summary ===
print(head(result_minimal$summary, 3))
#>            task_id  status rmse__value rmse__n_obs bias__value coverage__mean
#> 1 d001_f001_r00001 success    2.160240         100  0.42322202              1
#> 2 d001_f001_r00002 success    2.107857         100  0.42394617              1
#> 3 d001_f001_r00003 success    2.401193         100  0.09284687              1
#>   coverage__by_param__intercept coverage__by_param__slope
#> 1                             1                         1
#> 2                             1                         1
#> 3                             1                         1
#>   coverage__by_param__sigma timing_total rep_idx data_n data_intercept
#> 1                         1  0.033174515       1    100              1
#> 2                         1  0.049414396       2    100              1
#> 3                         1  0.005820751       3    100              1
#>   data_slope data_sigma fit_model
#> 1          2          1    linear
#> 2          2          1    linear
#> 3          2          1    linear

cat("\n=== Standard Retention Summary ===\n")
#> 
#> === Standard Retention Summary ===
print(head(result_standard$summary, 3))
#>            task_id  status rmse__value rmse__n_obs bias__value coverage__mean
#> 1 d001_f001_r00001 success    2.160240         100  0.42322202              1
#> 2 d001_f001_r00002 success    2.107857         100  0.42394617              1
#> 3 d001_f001_r00003 success    2.401193         100  0.09284687              1
#>   coverage__by_param__intercept coverage__by_param__slope
#> 1                             1                         1
#> 2                             1                         1
#> 3                             1                         1
#>   coverage__by_param__sigma rhat_max ess_bulk ess_tail divergent max_treedepth
#> 1                         1     1.01      400      350         0             0
#> 2                         1     1.01      400      350         0             0
#> 3                         1     1.01      400      350         0             0
#>   timing_total rep_idx data_n data_intercept data_slope data_sigma fit_model
#> 1  0.005964279       1    100              1          2          1    linear
#> 2  0.005661011       2    100              1          2          1    linear
#> 3  0.005592585       3    100              1          2          1    linear
```

Both profiles produce the same summary tibble with metric values. The
difference is in what additional data is available in the task results:

``` r
# Check what's available in task results (if any)
cat("Minimal retention - task results available:", 
    !is.null(result_minimal$task_results), "\n")
#> Minimal retention - task results available: TRUE
cat("Standard retention - task results available:", 
    !is.null(result_standard$task_results), "\n")
#> Standard retention - task results available: TRUE
```

### Memory Estimation

Use [`object.size()`](https://rdrr.io/r/utils/object.size.html) and
[`estimate_size()`](https://sims1253.github.io/bayesim/reference/estimate_size.md)
to check memory usage:

``` r
# Estimate size of the result objects
cat("Minimal result size:", 
    format(object.size(result_minimal), units = "KB"), "\n")
#> Minimal result size: 70.8 Kb
cat("Standard result size:", 
    format(object.size(result_standard), units = "KB"), "\n")
#> Standard result size: 89.8 Kb

# Estimate size of individual components
if (!is.null(result_standard$task_results) && length(result_standard$task_results) > 0) {
  first_task <- result_standard$task_results[[1]]
  cat("\nSingle task result size:", 
      format(estimate_size(first_task), units = "KB"), "\n")
}
#> 
#> Single task result size: 3752
```

## Tips for Memory Management

### Monitor Memory with `gc()` and `object.size()`

Regularly check memory usage during development:

``` r
# Before simulation
gc()
start_mem <- sum(gc()[, 2])  # Get used memory in MB

# Run simulation
result <- run_simulation(config, progress = FALSE)

# After simulation
gc()
end_mem <- sum(gc()[, 2])
cat("Memory used:", end_mem - start_mem, "MB\n")

# Check result size
result_size <- object.size(result)
cat("Result object size:", format(result_size, units = "MB"), "\n")
```

### Use Test Runs to Estimate Memory Needs

Before running a large study, do a small test run to estimate per-task
memory:

``` r
# Test run with 5 tasks
test_config <- simulation_config(
  data_grid = data.frame(n = 100, sigma = 1),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_data_generator,
  fitter = my_fitter,
  metrics = list(rmse_metric()),
  n_replicates = 5L,
  seed = 42L,
  retain = resolve_retention("debug")  # Test with full retention
)

test_result <- run_simulation(test_config, progress = FALSE)

# Estimate per-task memory
if (!is.null(test_result$task_results) && length(test_result$task_results) > 0) {
  per_task_size <- estimate_size(test_result$task_results[[1]])
  total_tasks <- 10000  # Your planned study size
  estimated_total <- per_task_size * total_tasks
  
  cat("Per-task size:", format(per_task_size, units = "KB"), "\n")
  cat("Estimated total for", total_tasks, "tasks:", 
      format(estimated_total, units = "GB"), "\n")
}
```

### Consider Splitting Very Large Studies

If your study is too large for available memory even with minimal
retention, split it into chunks:

``` r
# Split by data configuration
n_values <- c(100, 500, 1000)

for (n_val in n_values) {
  chunk_config <- simulation_config(
    data_grid = data.frame(n = n_val, sigma = c(0.5, 1, 2)),
    fit_grid = data.frame(model = c("m1", "m2")),
    data_generator = my_data_generator,
    fitter = my_fitter,
    metrics = list(rmse_metric()),
    n_replicates = 100L,
    seed = 42L,
    retain = resolve_retention("minimal"),
    result_path = paste0("chunk_n", n_val),
    checkpoint_every = 50L
  )
  
  result <- run_simulation(chunk_config, resume = "auto")
  
  # Save chunk summary
  saveRDS(result$summary, paste0("summary_n", n_val, ".rds"))
  
  # Explicit cleanup
  rm(result)
  gc()
}

# Combine summaries later
all_summaries <- lapply(n_values, function(n) {
  readRDS(paste0("summary_n", n, ".rds"))
})
combined_summary <- do.call(rbind, all_summaries)
```

### Use Checkpointing Strategically

Set `checkpoint_every` based on your model fitting time and available
memory:

``` r
# Fast models: checkpoint less frequently
fast_config <- simulation_config(
  # ... other params ...
  checkpoint_every = 500L  # Checkpoint every 500 tasks
)

# Slow models: checkpoint more frequently to avoid losing work
slow_config <- simulation_config(
  # ... other params ...
  checkpoint_every = 10L  # Checkpoint every 10 tasks
)
```

### Clean Up Old Checkpoints

After a simulation completes, you may want to remove old checkpoints to
free disk space:

``` r
# List checkpoints
checkpoint_ids <- list_checkpoints("my_simulation_results")
cat("Available checkpoints:", checkpoint_ids, "\n")

# Clean old checkpoints, keeping only the 3 most recent
clean_old_checkpoints("my_simulation_results", keep_n = 3)
```

## Summary

This vignette covered bayesim’s memory management system:

1.  **Retention profiles** control which artifacts are kept:

    - `"minimal"`: Only metrics (smallest memory)
    - `"standard"`: Metrics + diagnostics + warnings (default, balanced)
    - `"debug"`: Everything (largest memory, for debugging)

2.  **Custom retention** allows fine-grained control with explicit
    character vectors

3.  **Study size guidelines**:

    - Small (\< 100 tasks): Use `"debug"` or `"standard"`
    - Medium (100-1000 tasks): Use `"standard"`
    - Large (\> 1000 tasks): Use `"minimal"` or custom

4.  **Memory-bounded checkpointing** keeps memory usage bounded:

    - Results are checkpointed to disk periodically
    - Heavy objects are cleared from memory after checkpointing
    - `chunk_size` parameter controls when memory is freed

5.  **Monitoring tools** like
    [`object.size()`](https://rdrr.io/r/utils/object.size.html),
    [`estimate_size()`](https://sims1253.github.io/bayesim/reference/estimate_size.md),
    and [`gc()`](https://rdrr.io/r/base/gc.html) help you understand and
    plan for memory needs

Recommendations: - Start with `"standard"` retention for development -
Switch to `"minimal"` for runs with many tasks - Enable checkpointing
for studies \> 1000 tasks - Use `chunk_size` to control memory bounds
independently from checkpoint frequency - Test with small runs to
estimate memory needs before scaling up

For more information on checkpointing and resume functionality, see
[`vignette("getting-started")`](https://sims1253.github.io/bayesim/articles/getting-started.md).
For creating custom metrics that work efficiently with retention
profiles, see
[`vignette("custom-fitters")`](https://sims1253.github.io/bayesim/articles/custom-fitters.md).
