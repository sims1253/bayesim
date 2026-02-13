# bayesim (development version)

# bayesim 1.0.0

This is a major rewrite of bayesim, introducing a modern simulation framework 
with deterministic reproducibility, checkpoint/resume capabilities, and 
memory-bounded execution.

## Breaking Changes

* Complete rewrite of the simulation execution path
* New S7-based interface for custom fitters and metrics
* Old `full_simulation()`, `dataset_sim()`, `fit_sim()` functions are deprecated

## New Features

### Extension System

* `SimulationConfig` S7 class for validated, immutable configuration
* `Fitter` S7 abstract class for custom model fitting backends
* `Metric` S7 abstract class for custom metrics
* `MockFitter` for testing custom fitters and metrics

### Error Handling

* Structured error conditions: `bayesim_config_error`, `bayesim_contract_error`, 
  `bayesim_checkpoint_error`, `bayesim_internal_error` (fatal)
* Task-level recoverable errors: `bayesim_data_error`, `bayesim_fit_error`, 
  `bayesim_metric_error`
* Helper predicates: `is_bayesim_error()`, `is_fatal_error()`, `is_recoverable_error()`

### Result Types

* `bayesim_fit_result` S3 class for fit outputs
* `bayesim_task_result` S3 class for task outcomes  
* `bayesim_simulation_result` S3 class for complete runs

### Validation

* `validate_data_bundle()` for data generator output validation
* `validate_fitter_interface()` for fitter contract validation
* `validate_metric_interface()` for metric contract validation
* `validate_simulation_config()` for complete configuration validation

### Utilities

* Atomic file operations: `write_json_atomic()`, `write_rds_atomic()`
* Config fingerprinting: `compute_config_fingerprint()`
* Task ID formatting: `format_task_id()`, `parse_task_id()`
* Timing utilities: `make_timer()`
* Error capture: `capture_error_info()`

## Dependencies

* Added S7 for OOP
* Added digest for hashing
* Added jsonlite for serialization
* Added tibble for result data frames
* Moved brms, rstan, cmdstanr to Suggests

## Phase 2: Deterministic Engine Core

### Execution Engine

* `run_simulation()` main entry point for simulation runs
* `create_task_grid()` generates deterministic task table with precomputed RNG streams
* Task IDs in format `dXXX_fXXX_rXXXXX` for lexicographic ordering
* `execute_tasks()` iterates through tasks with progress bar and error tracking

### RNG Management

* `setup_global_rng()` initializes L'Ecuyer-CMRG RNG
* `set_task_rng()` restores per-task RNG state deterministically
* `advance_rng_stream()` pure function for advancing RNG state
* Each task gets independent, precomputed RNG stream

### Worker Execution

* `run_task()` executes single task: data generation, fitting, metrics
* `run_task_safe()` wrapper with fatal error propagation
* `build_metric_context()` precomputes context (predictions, log_lik, loo)
* `compute_all_metrics()` with required vs optional metric handling
* `apply_retention()` removes large objects based on retention policy

### Metric Registry

* `register_metric()` adds metrics with Metric subtype enforcement
* `get_metric()` retrieves metrics by name
* `list_metrics()` returns all registered metric names
* `unregister_metric()` removes metrics
* `clear_registry()` for testing (internal)

### Bug Fixes

* Fixed `set_task_rng()` to use explicit environment assignment
* Fixed `execute_tasks()` to pass proper S7 config object
* Fixed `run_task_safe()` to propagate fatal errors
* Removed duplicate `create_task_rng_streams()` function

## Phase 3: Checkpoint/Resume

### Checkpoint System

* `init_checkpoint_dir()` creates checkpoint directory structure
* `write_checkpoint()` atomically writes checkpoint with validation
* `read_checkpoint()` reads checkpoint with checksum verification
* `list_checkpoints()` lists available checkpoint IDs
* `get_latest_valid_checkpoint()` finds newest valid checkpoint
* Schema versioning for format compatibility
* Read-back validation for integrity

### Resume Logic

* `can_resume()` checks for valid resumable run
* `load_for_resume()` loads previous state with validation
* `get_resume_summary()` shows resumption summary
* `merge_task_grid_status()` merges task status from checkpoint
* `merge_results()` deduplicates results by task_id

### Integration

* `run_simulation()` supports `resume` and `force_restart` parameters
* Periodic checkpointing during execution
* Resume continues from pending tasks
* Prior results carried into final output

### Bug Fixes

* Fixed checkpoint write to use atomic RDS operations
* Fixed get_next_checkpoint_id() to ignore .tmp directories
* Fixed resume to validate both schema versions
* Fixed resume to not reinitialize manifest
* Fixed resume to carry prior results into final checkpoint

## Phase 4: Memory Management

### Retention System

* `resolve_retention()` resolves retention profiles to explicit options
* `apply_fit_retention()` removes non-retained fields from fit results
* `apply_task_retention()` adds retained fields to task results
* `estimate_size()` estimates object memory size
* `exceeds_size_threshold()` checks if result exceeds size limit
* `externalize_artifact()` moves large artifacts to external files

### Retention Profiles

* `minimal`: metrics only
* `standard`: metrics + diagnostics (default)
* `debug`: all fields retained

## Phase 5: brms Integration and Documentation

### BrmsFitter

* `BrmsFitter` class extending Fitter
* Supports rstan and cmdstanr backends
* Configurable MCMC settings (chains, iter, warmup, etc.)
* Full method implementations: fit, extract_draws, predict, log_lik, loo, diagnostics
* Automatic warning capture and diagnostic extraction

### Built-in Metrics

* `RmseMetric`: Root Mean Square Error
* `BiasMetric`: Mean prediction bias
* `CoverageMetric`: Credible interval coverage
* `PosteriorMeanMetric`: Posterior mean estimates
* `register_built_in_metrics()` for auto-registration

### Documentation

* Getting started vignette with examples
* Collate field added for proper file load order
