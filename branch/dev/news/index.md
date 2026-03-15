# Changelog

## bayesim 1.0.1

This release finalizes the rewrite follow-up work around the new
simulation API, resume/checkpoint behavior, and package documentation.

### Changes

- Removed legacy code that was superseded by the 1.0 rewrite:
  `simulation.R`, `inverse_forward_sampling.R`, `loo_handler.R`,
  `metric_list_handler.R`, `metric_lookup.R`, `ifs_sbc.R`,
  `ll_lookup.R`, `prefit.R`, `parallel_helpers.R`, and
  `simulation_building_blocks.R`.

- Removed corresponding man pages for all deleted functions.

- Simplified contracts, checkpoint, retention, worker, and simulation
  config code to eliminate dead paths and tighten validation.

- Updated NAMESPACE and DESCRIPTION to reflect the reduced API surface.

- Code quality improvements from desloppify review.

- Standardized the public workflow around
  [`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md),
  [`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md),
  and
  [`resume_simulation()`](https://sims1253.github.io/bayesim/reference/resume_simulation.md).

- Added support for explicit `task_grid`, `chunk_size`, conditional
  retention, manifest-based resume, stricter duplicate-check handling,
  and future-based batch execution.

- Switched the default
  [`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
  backend to `"cmdstanr"`.

- Kept `checkpoint_format = "rds"` as the supported checkpoint backend
  and made unsupported `"parquet"` requests fail fast.

- Externalized large metric payloads into artifacts to avoid excessively
  wide summary tables.

- Updated README, vignettes, roxygen docs, man pages, and NAMESPACE to
  match the rewritten API and current package behavior.

## bayesim 1.0.0

This is a major rewrite of bayesim, introducing a modern simulation
framework with deterministic reproducibility, checkpoint/resume
capabilities, and memory-bounded execution.

### Breaking Changes

- Complete rewrite of the simulation execution path
- New S7-based interface for custom fitters and metrics
- Old
  [`full_simulation()`](https://sims1253.github.io/bayesim/reference/full_simulation.md),
  [`dataset_sim()`](https://sims1253.github.io/bayesim/reference/dataset_sim.md),
  [`fit_sim()`](https://sims1253.github.io/bayesim/reference/fit_sim.md)
  functions are deprecated

### New Features

#### Extension System

- `SimulationConfig` S7 class for validated, immutable configuration
- `Fitter` S7 abstract class for custom model fitting backends
- `Metric` S7 abstract class for custom metrics
- `MockFitter` for testing custom fitters and metrics

#### Error Handling

- Structured error conditions: `bayesim_config_error`,
  `bayesim_contract_error`, `bayesim_checkpoint_error`,
  `bayesim_internal_error` (fatal)
- Task-level recoverable errors: `bayesim_data_error`,
  `bayesim_fit_error`, `bayesim_metric_error`
- Helper predicates:
  [`is_bayesim_error()`](https://sims1253.github.io/bayesim/reference/is_bayesim_error.md),
  [`is_fatal_error()`](https://sims1253.github.io/bayesim/reference/is_fatal_error.md),
  [`is_recoverable_error()`](https://sims1253.github.io/bayesim/reference/is_recoverable_error.md)

#### Result Types

- `bayesim_fit_result` S3 class for fit outputs
- `bayesim_task_result` S3 class for task outcomes  
- `bayesim_simulation_result` S3 class for complete runs

#### Validation

- [`validate_data_bundle()`](https://sims1253.github.io/bayesim/reference/validate_data_bundle.md)
  for data generator output validation
- [`validate_fitter_interface()`](https://sims1253.github.io/bayesim/reference/validate_fitter_interface.md)
  for fitter contract validation
- [`validate_metric_interface()`](https://sims1253.github.io/bayesim/reference/validate_metric_interface.md)
  for metric contract validation
- [`validate_simulation_config()`](https://sims1253.github.io/bayesim/reference/validate_simulation_config.md)
  for complete configuration validation

#### Utilities

- Atomic file operations:
  [`write_json_atomic()`](https://sims1253.github.io/bayesim/reference/write_json_atomic.md),
  [`write_rds_atomic()`](https://sims1253.github.io/bayesim/reference/write_rds_atomic.md)
- Config fingerprinting:
  [`compute_config_fingerprint()`](https://sims1253.github.io/bayesim/reference/compute_config_fingerprint.md)
- Task ID formatting:
  [`format_task_id()`](https://sims1253.github.io/bayesim/reference/format_task_id.md),
  [`parse_task_id()`](https://sims1253.github.io/bayesim/reference/parse_task_id.md)
- Timing utilities:
  [`make_timer()`](https://sims1253.github.io/bayesim/reference/make_timer.md)
- Error capture:
  [`capture_error_info()`](https://sims1253.github.io/bayesim/reference/capture_error_info.md)

### Dependencies

- Added S7 for OOP
- Added digest for hashing
- Added jsonlite for serialization
- Added tibble for result data frames
- Moved brms, rstan, cmdstanr to Suggests

### Phase 2: Deterministic Engine Core

#### Execution Engine

- [`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
  main entry point for simulation runs
- [`create_task_grid()`](https://sims1253.github.io/bayesim/reference/create_task_grid.md)
  generates deterministic task table with precomputed RNG streams
- Task IDs in format `dXXX_fXXX_rXXXXX` for lexicographic ordering
- [`execute_tasks()`](https://sims1253.github.io/bayesim/reference/execute_tasks.md)
  iterates through tasks with progress bar and error tracking

#### RNG Management

- [`setup_global_rng()`](https://sims1253.github.io/bayesim/reference/setup_global_rng.md)
  initializes L’Ecuyer-CMRG RNG
- [`set_task_rng()`](https://sims1253.github.io/bayesim/reference/set_task_rng.md)
  restores per-task RNG state deterministically
- [`advance_rng_stream()`](https://sims1253.github.io/bayesim/reference/advance_rng_stream.md)
  pure function for advancing RNG state
- Each task gets independent, precomputed RNG stream

#### Worker Execution

- [`run_task()`](https://sims1253.github.io/bayesim/reference/run_task.md)
  executes single task: data generation, fitting, metrics
- [`run_task_safe()`](https://sims1253.github.io/bayesim/reference/run_task_safe.md)
  wrapper with fatal error propagation
- [`build_metric_context()`](https://sims1253.github.io/bayesim/reference/build_metric_context.md)
  precomputes context (predictions, log_lik, loo)
- [`compute_all_metrics()`](https://sims1253.github.io/bayesim/reference/compute_all_metrics.md)
  with required vs optional metric handling
- [`apply_retention()`](https://sims1253.github.io/bayesim/reference/apply_retention.md)
  removes large objects based on retention policy

#### Metric Registry

- `register_metric()` adds metrics with Metric subtype enforcement
- `get_metric()` retrieves metrics by name
- `list_metrics()` returns all registered metric names
- `unregister_metric()` removes metrics
- [`clear_registry()`](https://sims1253.github.io/bayesim/reference/clear_registry.md)
  for testing (internal)

#### Bug Fixes

- Fixed
  [`set_task_rng()`](https://sims1253.github.io/bayesim/reference/set_task_rng.md)
  to use explicit environment assignment
- Fixed
  [`execute_tasks()`](https://sims1253.github.io/bayesim/reference/execute_tasks.md)
  to pass proper S7 config object
- Fixed
  [`run_task_safe()`](https://sims1253.github.io/bayesim/reference/run_task_safe.md)
  to propagate fatal errors
- Removed duplicate
  [`create_task_rng_streams()`](https://sims1253.github.io/bayesim/reference/create_task_rng_streams.md)
  function

### Phase 3: Checkpoint/Resume

#### Checkpoint System

- [`init_checkpoint_dir()`](https://sims1253.github.io/bayesim/reference/init_checkpoint_dir.md)
  creates checkpoint directory structure
- [`write_checkpoint()`](https://sims1253.github.io/bayesim/reference/write_checkpoint.md)
  atomically writes checkpoint with validation
- [`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md)
  reads checkpoint with checksum verification
- [`list_checkpoints()`](https://sims1253.github.io/bayesim/reference/list_checkpoints.md)
  lists available checkpoint IDs
- [`get_latest_valid_checkpoint()`](https://sims1253.github.io/bayesim/reference/get_latest_valid_checkpoint.md)
  finds newest valid checkpoint
- Schema versioning for format compatibility
- Read-back validation for integrity

#### Resume Logic

- [`can_resume()`](https://sims1253.github.io/bayesim/reference/can_resume.md)
  checks for valid resumable run
- [`load_for_resume()`](https://sims1253.github.io/bayesim/reference/load_for_resume.md)
  loads previous state with validation
- [`get_resume_summary()`](https://sims1253.github.io/bayesim/reference/get_resume_summary.md)
  shows resumption summary
- [`merge_task_grid_status()`](https://sims1253.github.io/bayesim/reference/merge_task_grid_status.md)
  merges task status from checkpoint
- [`merge_results()`](https://sims1253.github.io/bayesim/reference/merge_results.md)
  deduplicates results by task_id

#### Integration

- [`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
  supports `resume` and `force_restart` parameters
- Periodic checkpointing during execution
- Resume continues from pending tasks
- Prior results carried into final output

#### Bug Fixes

- Fixed checkpoint write to use atomic RDS operations
- Fixed get_next_checkpoint_id() to ignore .tmp directories
- Fixed resume to validate both schema versions
- Fixed resume to not reinitialize manifest
- Fixed resume to carry prior results into final checkpoint

### Phase 4: Memory Management

#### Retention System

- [`resolve_retention()`](https://sims1253.github.io/bayesim/reference/resolve_retention.md)
  resolves retention profiles to explicit options
- [`apply_fit_retention()`](https://sims1253.github.io/bayesim/reference/apply_fit_retention.md)
  removes non-retained fields from fit results
- [`apply_task_retention()`](https://sims1253.github.io/bayesim/reference/apply_task_retention.md)
  adds retained fields to task results
- [`estimate_size()`](https://sims1253.github.io/bayesim/reference/estimate_size.md)
  estimates object memory size
- [`exceeds_size_threshold()`](https://sims1253.github.io/bayesim/reference/exceeds_size_threshold.md)
  checks if result exceeds size limit
- [`externalize_artifact()`](https://sims1253.github.io/bayesim/reference/externalize_artifact.md)
  moves large artifacts to external files

#### Retention Profiles

- `minimal`: metrics only
- `standard`: metrics + diagnostics (default)
- `debug`: all fields retained

### Phase 5: brms Integration and Documentation

#### BrmsFitter

- `BrmsFitter` class extending Fitter
- Supports rstan and cmdstanr backends
- Configurable MCMC settings (chains, iter, warmup, etc.)
- Full method implementations: fit, extract_draws, predict, log_lik,
  loo, diagnostics
- Automatic warning capture and diagnostic extraction

#### Built-in Metrics

- `RmseMetric`: Root Mean Square Error
- `BiasMetric`: Mean prediction bias
- `CoverageMetric`: Credible interval coverage
- `PosteriorMeanMetric`: Posterior mean estimates
- [`register_built_in_metrics()`](https://sims1253.github.io/bayesim/reference/register_built_in_metrics.md)
  for auto-registration

#### Documentation

- Getting started vignette with examples
- Collate field added for proper file load order
