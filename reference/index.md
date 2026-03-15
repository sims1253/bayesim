# Package index

## Configuration

Functions for configuring simulation studies

- [`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
  : Create a Simulation Configuration
- [`is_simulation_config()`](https://sims1253.github.io/bayesim/reference/is_simulation_config.md)
  : Check if Object is a SimulationConfig
- [`validate_config_completeness()`](https://sims1253.github.io/bayesim/reference/validate_config_completeness.md)
  : Validate SimulationConfig Completeness
- [`get_total_tasks()`](https://sims1253.github.io/bayesim/reference/get_total_tasks.md)
  : Get Total Number of Tasks in Configuration
- [`compute_config_fingerprint()`](https://sims1253.github.io/bayesim/reference/compute_config_fingerprint.md)
  : Compute Configuration Fingerprint
- [`as_config_spec()`](https://sims1253.github.io/bayesim/reference/as_config_spec.md)
  : Convert SimulationConfig to Plain List for Hashing

## Running Simulations

Functions for executing simulation studies

- [`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
  : Run a Simulation Study
- [`resume_simulation()`](https://sims1253.github.io/bayesim/reference/resume_simulation.md)
  : Resume a simulation from an existing result directory
- [`can_resume()`](https://sims1253.github.io/bayesim/reference/can_resume.md)
  : Check if Resumable Run Exists

## Fitters

Classes for fitting Bayesian models

- [`Fitter()`](https://sims1253.github.io/bayesim/reference/Fitter.md) :
  Abstract Fitter Class
- [`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
  : Brms Fitter
- [`MockFitter()`](https://sims1253.github.io/bayesim/reference/MockFitter.md)
  : Mock Fitter for Testing (Internal)
- [`fit()`](https://sims1253.github.io/bayesim/reference/fit.md) : Fit a
  Bayesian Model
- [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)
  : Extract Posterior Draws
- [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
  : Generate Predictions
- [`log_lik()`](https://sims1253.github.io/bayesim/reference/log_lik.md)
  : Compute Pointwise Log-Likelihood
- [`loo()`](https://sims1253.github.io/bayesim/reference/loo.md) :
  Compute LOO-CV
- [`diagnostics()`](https://sims1253.github.io/bayesim/reference/diagnostics.md)
  : Extract Fit Diagnostics
- [`validate_fitter()`](https://sims1253.github.io/bayesim/reference/validate_fitter.md)
  : Validate a Fitter Object

## Metrics

Classes and built-in metrics for computing simulation results

- [`Metric()`](https://sims1253.github.io/bayesim/reference/Metric.md) :
  Metric Abstract Class
- [`compute()`](https://sims1253.github.io/bayesim/reference/compute.md)
  : Compute Metric Values
- [`RmseMetric()`](https://sims1253.github.io/bayesim/reference/RmseMetric.md)
  [`rmse_metric()`](https://sims1253.github.io/bayesim/reference/RmseMetric.md)
  : RMSE Metric
- [`BiasMetric()`](https://sims1253.github.io/bayesim/reference/BiasMetric.md)
  [`bias_metric()`](https://sims1253.github.io/bayesim/reference/BiasMetric.md)
  : Bias Metric
- [`CoverageMetric()`](https://sims1253.github.io/bayesim/reference/CoverageMetric.md)
  [`coverage_metric()`](https://sims1253.github.io/bayesim/reference/CoverageMetric.md)
  : Coverage Metric
- [`PosteriorMeanMetric()`](https://sims1253.github.io/bayesim/reference/PosteriorMeanMetric.md)
  [`posterior_mean_metric()`](https://sims1253.github.io/bayesim/reference/PosteriorMeanMetric.md)
  : Posterior Mean Metric
- [`validate_metric_output()`](https://sims1253.github.io/bayesim/reference/validate_metric_output.md)
  : Validate Metric Output
- [`flatten_metric_output()`](https://sims1253.github.io/bayesim/reference/flatten_metric_output.md)
  : Flatten Metric Output

## Results

Functions for working with simulation results

- [`is_bayesim_fit_result()`](https://sims1253.github.io/bayesim/reference/is_bayesim_fit_result.md)
  : Check if object is a bayesim_fit_result
- [`validate_bayesim_fit_result()`](https://sims1253.github.io/bayesim/reference/validate_bayesim_fit_result.md)
  : Validate a bayesim_fit_result object
- [`new_fit_result()`](https://sims1253.github.io/bayesim/reference/new_fit_result.md)
  : Create a new bayesim_fit_result object
- [`new_task_result()`](https://sims1253.github.io/bayesim/reference/new_task_result.md)
  : Create a new bayesim_task_result object
- [`new_simulation_result()`](https://sims1253.github.io/bayesim/reference/new_simulation_result.md)
  : Create a new bayesim_simulation_result object

## Checkpointing

Functions for checkpoint and resume functionality

- [`init_checkpoint_dir()`](https://sims1253.github.io/bayesim/reference/init_checkpoint_dir.md)
  : Initialize checkpoint directory
- [`write_checkpoint()`](https://sims1253.github.io/bayesim/reference/write_checkpoint.md)
  : Write checkpoint
- [`read_checkpoint()`](https://sims1253.github.io/bayesim/reference/read_checkpoint.md)
  : Read checkpoint
- [`get_resume_summary()`](https://sims1253.github.io/bayesim/reference/get_resume_summary.md)
  : Get Resume Summary
- [`load_for_resume()`](https://sims1253.github.io/bayesim/reference/load_for_resume.md)
  : Load Run for Resume

## Retention

Control what artifacts are retained in results

- [`resolve_retention()`](https://sims1253.github.io/bayesim/reference/resolve_retention.md)
  : Resolve retention specification
- [`RETENTION_PROFILES`](https://sims1253.github.io/bayesim/reference/RETENTION_PROFILES.md)
  : Retention profiles

## Utilities

Helper functions and utilities

- [`` `%||%` ``](https://sims1253.github.io/bayesim/reference/null-coalescing.md)
  : Null coalescing operator
- [`cluster_setup()`](https://sims1253.github.io/bayesim/reference/cluster_setup.md)
  : Convenience Function to set up a cluster used for multiprocessing
- [`ifs_SBC()`](https://sims1253.github.io/bayesim/reference/ifs_SBC.md)
  : Full inverse forward sampling supported SBC
