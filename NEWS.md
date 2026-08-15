# bayesim (development version)

Post-review hardening of the 2.0.0 engine, metrics, and analysis layer.

## Engine and resume

* Fixed the legacy-resume truth/diagnostics round-trip: resumed runs no
  longer lose or mangle recorded truths and fit diagnostics when prior task
  results are reloaded from a checkpoint.
* Fatal mid-batch errors now persist the batch's already-completed sibling
  outcomes before re-raising, so a crash no longer discards finished work.
* Resume no longer double-loads the full run history (prior results were
  loaded once to resume and again during the run).
* In-memory (`result_path = NULL`) run-store writes are now linear in the
  number of completed tasks instead of repeatedly rewriting the full state.
* Adaptive stopping now warns when its evaluation step fails instead of
  passing silently.

## Metrics

* Removed the unvalidated mori shared-memory model-bank integration; model
  banks travel to daemons by ordinary serialization.
* Removed `rstar_metric()` (and with it the caret/randomForest
  dependencies) and the `rmse_test_metric()` alias, and merged
  `convergence_metric()` into an extended `sampler_diagnostics_metric()`
  (now emitting `rhat_max`, `ess_bulk_min`, `ess_tail_min`, `divergent`,
  and `max_treedepth`).
* `pos_prob_metric()`'s `by_param` field now declares mean/sd aggregation
  instead of a binomial MCSE, which was wrong for a posterior probability
  mean.
* Metric NA-degradation paths now emit schema-conformant fields (present
  with `NA` values) instead of dropping fields from the flattened summary.

## Fitters and errors

* New `supports_epred` fitter capability gating `predict_epred()` — `TRUE`
  for `LinearRegressionFitter()`/`BrmsFitter()`, dynamic (set when an
  `epred` generated quantity is declared) for `CmdStanFitter()`.
* `bayesim_contract_error()` is now exported.

## Analysis and reporting

* `report()` was renamed to `render_report()` to stop colliding with the
  generic of the easystats *report* package. `report()` remains as a
  deprecated alias that forwards to `render_report()` and warns once per
  session.

# bayesim 2.0.0

A ground-up rewrite of the simulation engine, fitters, generators, metrics, and
analysis layer, redesigned around the needs of simulation-method studies
(Morris, White & Crowther 2019). Breaking across the public API; the package
remains GitHub-only and lifecycle-experimental.

## Review hardening

* Fixed default simulation summaries with failed tasks: error payloads and
  other flattened metric metadata no longer become condition columns, and
  each aggregate now reports the number of finite values used.
* Adaptive stopping now requires the requested MCSE target in every condition
  cell, not only the first.
* Retained predictions are computed and returned even without a prediction
  metric; checkpoint lightening now honors all explicit retention options and
  preserves truths.
* `run_simulation()` restores the caller's RNG state and kind. Stochastic
  metrics use stable metric-specific sub-seeds, so metric order is irrelevant.
* Task-grid bookkeeping is vectorized by batch. The run store is append-only:
  completed outcomes are written once as immutable, mirrored outcome shards,
  and each checkpoint commit directory records its deltas plus metadata.
  `keep_checkpoints = 2L` by default prunes old checkpoint commits, keeping
  one fallback commit; immutable outcome/ledger history is never pruned, so
  durable storage grows roughly linearly with completed tasks.
* Fixed S7 class checks that had made automatic BrmsFitter model-bank and
  generator shortcuts unreachable. Precompiled brms models now warn when
  explicit priors are missing, and preflight reports deduplicated compile
  counts.
* Added exported `config_fingerprint()`, repaired the targets/HPC/SBC
  vignettes, consolidated `rmse_test_metric()` onto `pred_rmse_metric()`, and
  removed unused dependencies, helpers, and a committed Stan binary.

## Fitters and the fitter contract

* **Contract matrices are now `S x N` (draws x observations) everywhere** —
  `log_lik_matrix`, `predict_fit()$predicted_samples`, `predict_epred`. The
  brms/loo convention, enforced by non-square `validate_fitter()` smoke tests.
* **New `LinearRegressionFitter()`** — an exact conjugate Normal-Inverse-Gamma
  Bayesian linear regression fitter. Real posteriors in milliseconds with zero
  Stan; the package's executable-docs teaching backbone.
* **New `CmdStanFitter()`** — run user-supplied Stan programs via cmdstanr
  (declared `log_lik`/`epred` generated quantities). No model-bank integration
  (cmdstanr caches binaries by file hash).
* **New `brms_model()` / `model_grid()`** — assemble tidy `fit_grid`s of brms
  model specs (formula/family/prior list-columns), validated at construction.
* `MockFitter` demoted to internal (testing only).
* `BrmsFitter` diagnostics now computed over **all** parameters (group-level,
  distributional, sigma) via `posterior::summarise_draws`, not just fixed
  effects; `loo_fit` uses chain-aware `r_eff`, matching `brms::loo()`.
* Consolidated fitter/metric validation into the exported `validate_fitter()`
  / `validate_metric()` (removed `check_fitter_class`,
  `validate_fitter_interface`, `validate_metric_interface`).

## Renamed generics (no aliases — pre-release API)

| old | new |
|---|---|
| `fit()` | `fit_model()` |
| `compute()` | `compute_metric()` |
| `log_lik()` | `log_lik_matrix()` |
| `diagnostics()` | `fit_diagnostics()` |

bayesim no longer exports `fit`/`compute`/`log_lik`/`diagnostics`, so it stops
masking `generics::fit`, `dplyr::compute`, and `brms::log_lik`.

## Transport: purrr + mirai

* Dispatch is now a single `purrr::map()` + `purrr::in_parallel()` code path;
  mirai remains the daemon engine. `run_task_safe()` is total (fatal conditions
  are captured and re-raised after the batch with their full class chain),
  removing the cross-boundary condition-restoration machinery.
* `run_simulation(config, workers = N)` sets up and tears down mirai daemons
  for the run (the simple path); `mirai::daemons()` remains for advanced/HPC.
* The model bank is now shared across local daemons via [mori](https://shikokuchuo.net/mori/)
  shared memory, so each daemon zero-copy maps the bank instead of
  deserializing a private copy. Locality is auto-detected from the mirai
  daemon URL (machine-local transports and loopback TCP qualify); remote and
  wildcard-bound daemons, and sequential runs, are passed through unchanged.
  `mori` is an optional suggested dependency because its current releases
  require R >= 4.3 while bayesim supports R >= 4.1. Without it, or when shared
  memory creation fails, bayesim safely uses ordinary serialized dispatch.

## Metrics and analysis (Morris et al. framing)

* **`performance_measures()`** — bias, empirical SE, MSE, coverage, average
  model SE, and `n_sim`, each with its Morris et al. MCSE, per estimand and
  condition. The centerpiece of the analysis layer.
* Task results now carry the data-generating truth, flattened to
  `truth__<param>` summary columns — enabling parameter-recovery analysis and
  `plot_recovery()` for generator-drawn truths.
* Prediction metrics renamed honestly: `rmse`/`bias`/`mae`/`mse` →
  `pred_*_metric`. They **refuse to silently fall back to the training set**:
  default to test data, NA with a one-time warning when no test set.
* `Metric` gained a `summary_type` property (`"mean"`, `"proportion"`,
  `"none"`), recorded on the result and consumed by `summarize_simulation()` —
  replacing both the 0-1-column coverage heuristic and the mathematically
  wrong per-task-RMSE delta-method MCSE.
* Metric context is computed on the **test set** when one exists: predictions
  and pointwise log-lik previously came from the training data while every
  consuming metric compared against the test response. LOO metrics
  (`rmse_loo`, `r2_loo`) now always compare against the training response
  (LOO is in-sample by construction).
* `true_params` / `vars_of_interest` are now optional (truth-free studies run).
* `rstar_metric()` reimplemented for real (was a placebo returning NA): uses
  `posterior::rstar()` on per-chain draws from the underlying fit; degrades to
  NA with a one-time warning for chain-less fitters.
* `plot_coverage()` redesigned as a point-range plot with MCSE error bars.

## Config knobs

* `chunk_size` and the deprecated `max_in_memory` merged into a single
  `checkpoint_every` knob.
* Data-generator signature is now `(data_spec, task_ctx)`; `task_ctx$seed`
  carries the integer for backends that need one.
* `retain`, `max_errors`, and `checkpoint_format` excluded from the config
  fingerprint — changing retention no longer invalidates resume.

## Runtime UX

* `preflight(config)` reports task count, grid shape, metrics' needs vs fitter
  capabilities, daemons status, compile count (auto one-liner in
  `run_simulation()`).
* `failed_tasks(result)` accessor; compact failure summary printed at run end.
* `print()` shows a metrics preview; `as_tibble.bayesim_simulation_result`
  returns the summary for tidyverse piping.
* Better `seed` error message.

## Error API

* Exported error constructors narrowed to `bayesim_data_error`,
  `bayesim_fit_error`, `bayesim_metric_error`, `bayesim_config_error` plus
  `is_bayesim_error`, `is_fatal_error`, `is_recoverable_error`. Base
  `bayesim_error`/contract/checkpoint/internal constructors are internal.

---

## Earlier 2.0.0 implementation milestones

### Engine
* Parallelization switched from `future`/`future.apply` to **mirai**. The
  engine owns parallelism; `SimulationConfig` gained a `daemon_setup` property
  for per-daemon initialization (e.g. cmdstan path configuration). Tasks run
  sequentially when no mirai daemons are set.
* Per-task RNG is L'Ecuyer-CMRG streams restored in the worker before the data
  generator runs; generators consume the ambient state (do not re-seed).

### BrmsFitter and the model bank
* `BrmsFitter` now compiles each distinct model spec **once** via
  `brms::brm(chains = 0)` at `run_simulation()` start and reuses the compiled
  binary via `stats::update(recompile = FALSE)` per task — eliminating the
  catastrophic per-task recompilation of 1.x.
* Added `precompile` (default `TRUE`) and `stan_args` properties. A structural
  data-mismatch guard aborts loudly if task data would require recompilation.
* `extract_brms_timings` now reports real warmup/sample timings from the
  stanfit object for both cmdstanr and rstan backends.

### Generators (new)
* `fixed_truth_generator()`, `prior_predictive_generator()`,
  `ifs_generator()` — factory constructors for the standard generator
  signature `(data_spec, task_ctx) -> data_bundle`. IFS and
  prior-predictive use a deterministic draw index (`task_ctx$rep_idx`) so SBC
  ranks are well-defined and resume is reproducible.
* Inverse forward sampling internals (`brms_full_ppred`,
  `brms_response_sequence`, `nodes_by_depth`) ported from 0.x and decoupled
  from SBC/bayeshear/future.

### Metrics (expanded)
* 14 new metric constructors: `mae_metric`, `mse_metric`, `pos_prob_metric`,
  `posterior_summary_metric`, `convergence_metric`, `sampler_diagnostics_metric`,
  `rank_metric`, `rstar_metric`, `elpd_loo_metric`, `rmse_loo_metric`,
  `r2_loo_metric`, `elpd_test_metric`, `rmse_test_metric`, `r2_test_metric`.
* Removed the internal metric registry; metrics are constructed directly.

### Analysis & reporting (new)
* `summarize_simulation()` — per-condition aggregation with Monte Carlo
  standard errors (rsimsum formulas).
* SBC diagnostics: `sbc_ranks()`, `plot_rank_hist()`, `plot_rank_ecdf()` (with
  simultaneous confidence bands), `plot_recovery()`, `plot_coverage()`,
  `plot_metric()`. Plotting requires ggplot2 (Suggests, loaded on demand).

## Breaking changes
* Public API contracted to a curated surface of 65 exports (see
  `_pkgdown.yml`).
* The `loo` S7 generic was renamed `loo_fit` to avoid clashing with
  `loo::loo()`. Custom fitters must implement `loo_fit`.
* Operators `%+%` and `%||%` are no longer exported (use `rlang::%||%` or the
  public equivalents).

## Reproducibility
Same seed + same package/backend versions + same platform produces identical
results. Across platforms (or backend version changes) results are
statistically equivalent but may differ in the least significant bits due to
floating-point and Stan version differences.

## Bug fixes (post-implementation review)
These fixes address defects found in a code review of the initial 2.0 build.
Each is verified by a behavioral test, not just `R CMD check` green.

* **IFS forward sampling (critical):** `ifs_generator()` no longer returns the
  pilot dataset unchanged. The 0.x response-dependency topological sort
  (`brms_response_sequence` → adjacency matrix → `nodes_by_depth`) is restored
  so the response column is actually simulated for univariate and multivariate
  models; `brms_full_ppred()` now takes a single draw and returns a single
  data.frame, eliminating the draw-index/list-index mismatch that produced
  `NULL` for `rep_idx > 1`.
* **vars_of_interest / draws-column mismatch (critical):** generators strip
  the `b_` prefix (`vars_of_interest = c("x","Intercept")`) but brms draws
  keep it (`c("b_x","b_Intercept","sigma")`), so every truth-comparing metric
  silently returned NA on real brms fits. A new `resolve_draw_columns()` helper
  maps cleaned names to actual columns (or errors) and is used by
  `coverage_metric`, `posterior_summary_metric`, `pos_prob_metric`,
  `posterior_summary_metric`, and `rank_metric`. Output field names stay on
  the cleaned names.
* **LOO-RMSE and LOO-R² (critical):** the invented formulas
  (`sqrt(-2*mean(elpd))`, `1-exp(2*elpd/n)`) are replaced with PSIS-based
  constructions on `loo::E_loo()`. `build_metric_context()` computes the PSIS
  object once (with chain-derived `r_eff`, falling back to `r_eff = NULL`);
  `rmse_loo` uses `E_loo(ppred, type="mean")` and reports max Pareto-k̂;
  `r2_loo` reproduces brms' `loo_R2()` variance construction verbatim and uses
  `posterior_epred` (a new `predict_epred` Fitter generic). Parity with
  `brms::loo_R2()` verified on a fixture fit.
* **SBC rank metric (critical for SBC validity):** `rank_metric()` gained
  `thin = "auto"` (default), thinning toward the min bulk-ESS across ranked
  variables to mitigate autocorrelation bias; `thin = FALSE` and integer strides
  are alternatives. Output now includes per-variable `n_ranks` (post-thinning
  sample size + 1); `sbc_ranks()` and `plot_rank_ecdf()` surface it.
* **BrmsFitter warning capture:** warnings are now captured in the `fit`
  method (covering both the fit and the lazy `summary()`-driven convergence
  warning from diagnostics extraction), for both the model-bank and fresh-compile
  paths. Previously the fresh path had no handler and warnings were lost.
* **mirai / model-bank efficiency:** the model bank and `daemon_setup` hook
  are shipped to daemons once per `run_simulation()` instead of per batch; the
  prefit-side Stan data-structure signature is computed once in
  `build_model_bank()` and cached (not recomputed per task); the session bank
  is cleared on run exit.
* **Cleanup:** `plot_rank_ecdf()` now uses true simultaneous ECDF bands
  (`adjust_gamma` ported from 0.x / Säilynoja et al. 2022) instead of an
  approximate KS bound; `posterior` moved from Suggests to Imports; the `dplyr`
  dependency was dropped (base-R row-binding replaces `bind_rows`);
  internal error constructors remain intentionally unexported; dead
  `hash_to_row` removed and a missing-`formula` guard added to
  `build_model_bank()`; the IFS bounds `resample` docstring now honestly
  describes the NA-out / truncate behavior and its rank-bias implication.
* **Metric flatten:** single-parameter named-vector metric outputs (e.g. a
  one-variable `rank`/`coverage` `by_param`) now carry the `__<param>` suffix
  when flattened, so downstream consumers grepping for
  `<metric>__by_param__<param>` find them.

# bayesim 1.0.1

This release finalizes the rewrite follow-up work around the new simulation API,
resume/checkpoint behavior, and package documentation.

## Changes

* Removed legacy code that was superseded by the 1.0 rewrite:
  `simulation.R`, `inverse_forward_sampling.R`, `loo_handler.R`,
  `metric_list_handler.R`, `metric_lookup.R`, `ifs_sbc.R`, `ll_lookup.R`,
  `prefit.R`, `parallel_helpers.R`, and `simulation_building_blocks.R`.
* Removed corresponding man pages for all deleted functions.
* Simplified contracts, checkpoint, retention, worker, and simulation config
  code to eliminate dead paths and tighten validation.
* Updated NAMESPACE and DESCRIPTION to reflect the reduced API surface.
* Code quality improvements from desloppify review.

* Standardized the public workflow around `simulation_config()`,
  `run_simulation()`, and `resume_simulation()`.
* Added support for explicit `task_grid`, `chunk_size`, conditional retention,
  manifest-based resume, stricter duplicate-check handling, and future-based
  batch execution.
* Switched the default `BrmsFitter()` backend to `"cmdstanr"`.
* Kept `checkpoint_format = "rds"` as the supported checkpoint backend and
  made unsupported `"parquet"` requests fail fast.
* Externalized large metric payloads into artifacts to avoid excessively wide
  summary tables.
* Updated README, vignettes, roxygen docs, man pages, and NAMESPACE to match
  the rewritten API and current package behavior.

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
