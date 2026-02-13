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
