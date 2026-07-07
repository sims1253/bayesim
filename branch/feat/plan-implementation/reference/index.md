# Package index

## Design and run studies

Configure and execute simulation studies

- [`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
  : Create a Simulation Configuration
- [`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
  : Run a Simulation Study
- [`resume_simulation()`](https://sims1253.github.io/bayesim/reference/resume_simulation.md)
  : Resume a simulation from an existing result directory
- [`preflight()`](https://sims1253.github.io/bayesim/reference/preflight.md)
  : Preflight check for a simulation configuration
- [`n_replicates_for_target()`](https://sims1253.github.io/bayesim/reference/n_replicates_for_target.md)
  : Required number of replicates for a target MCSE

## Fitters

Classes and generics for fitting Bayesian models

- [`Fitter()`](https://sims1253.github.io/bayesim/reference/Fitter.md) :
  Abstract Fitter Class
- [`LinearRegressionFitter()`](https://sims1253.github.io/bayesim/reference/LinearRegressionFitter.md)
  : Conjugate Bayesian Linear Regression Fitter
- [`BrmsFitter()`](https://sims1253.github.io/bayesim/reference/BrmsFitter.md)
  : Brms Fitter
- [`CmdStanFitter()`](https://sims1253.github.io/bayesim/reference/CmdStanFitter.md)
  : CmdStan Fitter (user-supplied Stan programs)
- [`brms_model()`](https://sims1253.github.io/bayesim/reference/brms_model.md)
  : Construct a single brms model specification
- [`model_grid()`](https://sims1253.github.io/bayesim/reference/model_grid.md)
  : Assemble a tibble of model specifications for a fit_grid
- [`fit_model()`](https://sims1253.github.io/bayesim/reference/fit_model.md)
  : Fit a Bayesian Model
- [`extract_draws()`](https://sims1253.github.io/bayesim/reference/extract_draws.md)
  : Extract Posterior Draws
- [`predict_fit()`](https://sims1253.github.io/bayesim/reference/predict_fit.md)
  : Generate Predictions
- [`predict_epred()`](https://sims1253.github.io/bayesim/reference/predict_epred.md)
  : Compute Posterior Expectation Predictions
- [`log_lik_matrix()`](https://sims1253.github.io/bayesim/reference/log_lik_matrix.md)
  : Compute Pointwise Log-Likelihood
- [`loo_fit()`](https://sims1253.github.io/bayesim/reference/loo_fit.md)
  : Compute LOO-CV
- [`fit_diagnostics()`](https://sims1253.github.io/bayesim/reference/fit_diagnostics.md)
  : Extract Fit Diagnostics

## Generators

Factory constructors for simulation data generators

- [`fixed_truth_generator()`](https://sims1253.github.io/bayesim/reference/fixed_truth_generator.md)
  : Construct a fixed-truth data generator
- [`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md)
  : Construct a prior-predictive data generator
- [`ifs_generator()`](https://sims1253.github.io/bayesim/reference/ifs_generator.md)
  : Construct an inverse forward sampling (IFS) data generator
- [`prior_draws_generator()`](https://sims1253.github.io/bayesim/reference/prior_draws_generator.md)
  : Construct a fitter-agnostic prior-draws data generator
- [`forward_sim_generator()`](https://sims1253.github.io/bayesim/reference/forward_sim_generator.md)
  : Construct a fitter-agnostic forward-simulation (IFS) data generator

## Metrics

Classes and built-in metrics for computing simulation results

- [`Metric()`](https://sims1253.github.io/bayesim/reference/Metric.md) :
  Metric Abstract Class
- [`compute_metric()`](https://sims1253.github.io/bayesim/reference/compute_metric.md)
  : Compute Metric Values
- [`PosteriorSummaryMetric()`](https://sims1253.github.io/bayesim/reference/PosteriorSummaryMetric.md)
  [`posterior_summary_metric()`](https://sims1253.github.io/bayesim/reference/PosteriorSummaryMetric.md)
  : Posterior Summary Metric
- [`CoverageMetric()`](https://sims1253.github.io/bayesim/reference/CoverageMetric.md)
  [`coverage_metric()`](https://sims1253.github.io/bayesim/reference/CoverageMetric.md)
  : Coverage Metric
- [`RankMetric()`](https://sims1253.github.io/bayesim/reference/RankMetric.md)
  [`rank_metric()`](https://sims1253.github.io/bayesim/reference/RankMetric.md)
  : SBC Rank Metric
- [`PosProbMetric()`](https://sims1253.github.io/bayesim/reference/PosProbMetric.md)
  [`pos_prob_metric()`](https://sims1253.github.io/bayesim/reference/PosProbMetric.md)
  : Positivity Probability Metric
- [`ConvergenceMetric()`](https://sims1253.github.io/bayesim/reference/ConvergenceMetric.md)
  [`convergence_metric()`](https://sims1253.github.io/bayesim/reference/ConvergenceMetric.md)
  : Convergence Metric
- [`SamplerDiagnosticsMetric()`](https://sims1253.github.io/bayesim/reference/SamplerDiagnosticsMetric.md)
  [`sampler_diagnostics_metric()`](https://sims1253.github.io/bayesim/reference/SamplerDiagnosticsMetric.md)
  : Sampler Diagnostics Metric
- [`RstarMetric()`](https://sims1253.github.io/bayesim/reference/RstarMetric.md)
  [`rstar_metric()`](https://sims1253.github.io/bayesim/reference/RstarMetric.md)
  : R\* Convergence Metric
- [`RmseMetric()`](https://sims1253.github.io/bayesim/reference/RmseMetric.md)
  [`pred_rmse_metric()`](https://sims1253.github.io/bayesim/reference/RmseMetric.md)
  : RMSE Metric
- [`BiasMetric()`](https://sims1253.github.io/bayesim/reference/BiasMetric.md)
  [`pred_bias_metric()`](https://sims1253.github.io/bayesim/reference/BiasMetric.md)
  : Bias Metric
- [`MaeMetric()`](https://sims1253.github.io/bayesim/reference/MaeMetric.md)
  [`pred_mae_metric()`](https://sims1253.github.io/bayesim/reference/MaeMetric.md)
  : MAE Metric
- [`MseMetric()`](https://sims1253.github.io/bayesim/reference/MseMetric.md)
  [`pred_mse_metric()`](https://sims1253.github.io/bayesim/reference/MseMetric.md)
  : MSE Metric
- [`ElpdLooMetric()`](https://sims1253.github.io/bayesim/reference/ElpdLooMetric.md)
  [`elpd_loo_metric()`](https://sims1253.github.io/bayesim/reference/ElpdLooMetric.md)
  : ELPD-LOO Metric
- [`RmseLooMetric()`](https://sims1253.github.io/bayesim/reference/RmseLooMetric.md)
  [`rmse_loo_metric()`](https://sims1253.github.io/bayesim/reference/RmseLooMetric.md)
  : RMSE-LOO Metric
- [`R2LooMetric()`](https://sims1253.github.io/bayesim/reference/R2LooMetric.md)
  [`r2_loo_metric()`](https://sims1253.github.io/bayesim/reference/R2LooMetric.md)
  : R-squared LOO Metric
- [`ElpdTestMetric()`](https://sims1253.github.io/bayesim/reference/ElpdTestMetric.md)
  [`elpd_test_metric()`](https://sims1253.github.io/bayesim/reference/ElpdTestMetric.md)
  : ELPD Test-Set Metric
- [`RmseTestMetric()`](https://sims1253.github.io/bayesim/reference/RmseTestMetric.md)
  [`rmse_test_metric()`](https://sims1253.github.io/bayesim/reference/RmseTestMetric.md)
  : RMSE Test-Set Metric
- [`R2TestMetric()`](https://sims1253.github.io/bayesim/reference/R2TestMetric.md)
  [`r2_test_metric()`](https://sims1253.github.io/bayesim/reference/R2TestMetric.md)
  : R-squared Test-Set Metric

## Analysis and reporting

Summaries, SBC diagnostics, and plotting functions for simulation
results

- [`summarize_simulation()`](https://sims1253.github.io/bayesim/reference/summarize_simulation.md)
  : Aggregate simulation results per condition
- [`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)
  : Simulation-method performance measures with Monte-Carlo standard
  errors
- [`failed_tasks()`](https://sims1253.github.io/bayesim/reference/failed_tasks.md)
  : Extract failed tasks from a simulation result
- [`sbc_ranks()`](https://sims1253.github.io/bayesim/reference/sbc_ranks.md)
  : Extract SBC ranks from a simulation result
- [`plot_rank_hist()`](https://sims1253.github.io/bayesim/reference/plot_rank_hist.md)
  : Plot SBC rank histograms
- [`plot_rank_ecdf()`](https://sims1253.github.io/bayesim/reference/plot_rank_ecdf.md)
  : Plot SBC rank ECDF with simultaneous confidence band
- [`plot_recovery()`](https://sims1253.github.io/bayesim/reference/plot_recovery.md)
  : Plot parameter recovery (truth vs posterior estimate)
- [`plot_coverage()`](https://sims1253.github.io/bayesim/reference/plot_coverage.md)
  : Plot coverage rates per condition/parameter
- [`plot_metric()`](https://sims1253.github.io/bayesim/reference/plot_metric.md)
  : Plot a metric across conditions
- [`report()`](https://sims1253.github.io/bayesim/reference/report.md) :
  Render a simulation-study report
- [`read_summary()`](https://sims1253.github.io/bayesim/reference/read_summary.md)
  : Read a simulation summary file

## Extension support

Validators, error constructors, and classifiers for extending bayesim

- [`validate_data_bundle()`](https://sims1253.github.io/bayesim/reference/validate_data_bundle.md)
  : Validate a Data Bundle
- [`validate_fitter()`](https://sims1253.github.io/bayesim/reference/validate_fitter.md)
  : Validate a Fitter Object
- [`validate_metric()`](https://sims1253.github.io/bayesim/reference/validate_metric.md)
  : Validate a Metric Object
- [`bayesim_config_error()`](https://sims1253.github.io/bayesim/reference/bayesim_config_error.md)
  : Configuration validation error (Fatal)
- [`bayesim_data_error()`](https://sims1253.github.io/bayesim/reference/bayesim_data_error.md)
  : Data generation/validation error (Recoverable)
- [`bayesim_fit_error()`](https://sims1253.github.io/bayesim/reference/bayesim_fit_error.md)
  : Model fitting error (Recoverable)
- [`bayesim_metric_error()`](https://sims1253.github.io/bayesim/reference/bayesim_metric_error.md)
  : Metric computation error (Recoverable)
- [`is_bayesim_error()`](https://sims1253.github.io/bayesim/reference/is_bayesim_error.md)
  : Check if a condition is a bayesim error
- [`is_fatal_error()`](https://sims1253.github.io/bayesim/reference/is_fatal_error.md)
  : Check if an error is fatal
- [`is_recoverable_error()`](https://sims1253.github.io/bayesim/reference/is_recoverable_error.md)
  : Check if an error is recoverable (task-level)
