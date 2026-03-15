# bayesim 1.0 Rewrite Plan (Revised)

**Version:** 2.1  
**Status:** Revised implementation blueprint  
**Target:** CRAN-ready package with deterministic, resumable, memory-bounded simulation engine

---

## Bottom Line

This rewrite keeps the architecture simple: **S7 only for extension contracts** (`SimulationConfig`, `Fitter`, `Metric`), **S3 for all runtime and persisted data**.

The core engine must be deterministic across sequential/parallel/resume modes by assigning each task an independent precomputed L'Ecuyer stream, propagating explicit seeds into backend samplers, and never depending on scheduler order.

Checkpointing should be append-safe, atomic on the same filesystem, and restartable using a validated ledger + chunked result files with strict schema versioning.

The plan must promise **reproducibility within documented backend/platform limits**, not unconditional bitwise equality across all Bayesian backends.

The authoritative on-disk checkpoint format must not depend on an optional package; columnar backends are an optimization layer, not a correctness dependency.

---

## Action Plan

1. Lock contracts first (interfaces, schemas, error classes, checkpoint/storage format, backend seed rules)
2. Build deterministic task execution + ledger before implementing advanced metrics
3. Add checkpoint/resume with corruption recovery, duplicate detection, and config fingerprint checks
4. Implement retention/memory controls and artifact externalization; benchmark optional columnar sinks against the core format
5. Finish with integration + stress tests (determinism, interruption, memory growth, storage tradeoffs)

**Effort estimate:** Medium (1-2d) for contract scaffolding; Large (1-2 weeks) for a full production implementation with robust testing, backend-specific hardening, and checkpoint recovery coverage.

---

## 1. Scope and Non-Goals

### In Scope
- Full rewrite of simulation execution path
- Formal extension contracts for custom fitters/metrics/data generators
- Deterministic reproducibility across:
  - `future::sequential`
  - `future::multisession`
  - interrupted/resumed runs
- File-based checkpoint and resume
- Memory-bounded execution with configurable artifact retention
- Test strategy including stress and failure-injection tests
- Backend-qualified reproducibility policy
- Storage backend benchmark and decision for checkpoint/result persistence

### Out of Scope (for v1 rewrite)
- Distributed cluster orchestration beyond `future`
- Backward compatibility with legacy object formats
- Automatic migration from old checkpoint formats

---

## 2. Architecture: Strict S7/S3 Split

### S7 (contracts and user-facing configuration)
Use S7 for:
- `SimulationConfig` (validated immutable-ish config object)
- `Fitter` (abstract extension point)
- `Metric` (abstract extension point)

### S3 (runtime state and persisted outputs)
Use S3 for:
- `bayesim_fit_result`
- `bayesim_task_result`
- `bayesim_simulation_result`
- ledger/result rows and checkpoint payloads

### Rule
- **S7 objects are never serialized directly into checkpoints**
- Before execution, S7 config is converted to a plain normalized list (`config_spec`) for hashing, persistence, and worker transport

---

## 3. Public API (Final)

```r
# Main entry point
run_simulation <- function(config, resume = c("auto", "never", "must"), progress = TRUE) {}

# Resume from an existing result directory
resume_simulation <- function(result_path, config = NULL, progress = TRUE) {}

# Config constructor
simulation_config <- function(
  data_grid = NULL,
  fit_grid = NULL,
  task_grid = NULL,
  data_generator,      # function(data_spec, seed, task_ctx) -> data_bundle
  fitter,              # S7 Fitter instance
  metrics,             # list of Metric objects
  n_replicates = 1L,
  seed,
  result_path = NULL,
  checkpoint_format = c("rds", "parquet"),
  checkpoint_every = 50L,
  chunk_size = NULL,
  retain = c("metrics", "diagnostics"),
  max_errors = Inf
) {}

# Built-in metric constructors
rmse_metric <- function(...) {}
bias_metric <- function(...) {}
coverage_metric <- function(...) {}
posterior_mean_metric <- function(...) {}
```

### API policy
- `run_simulation()` always starts from an explicit `SimulationConfig`; `resume = "auto"` resumes when compatible state already exists at `result_path`
- `resume_simulation()` is the ergonomic restart path after an interrupted session
- If the original run used only manifest-rehydratable components, `resume_simulation(result_path)` can resume from disk alone
- If the original run used custom closures or non-rehydratable objects, `resume_simulation()` must require `config` and fail with a clear error otherwise
- Built-in metrics are ordinary constructors, not names looked up from global package state
- `task_grid` is an escape hatch for sparse or prefiltered task plans; if supplied, it bypasses automatic Cartesian crossing of `data_grid x fit_grid x n_replicates`

### Defaults
- `resume = "auto"`
- `retain = c("metrics", "diagnostics")` as a shorthand for retaining the same artifacts for all task outcomes
- `checkpoint_format = "rds"` for the authoritative resume/checkpoint path; `"parquet"` is optional and only enabled when a supported columnar backend is available
- `checkpoint_every = 50L`
- `chunk_size = checkpoint_every` unless overridden
- `max_errors = Inf` (user can fail-fast by setting finite threshold)

### Input rules
- Exactly one of the following must be provided:
  - `task_grid`
  - both `data_grid` and `fit_grid`
- If `task_grid` is omitted, the engine constructs the task plan as the Cartesian product of `data_grid`, `fit_grid`, and `n_replicates`
- `metrics` must be a list of `Metric` objects; no exported global registry API
- `retain` may be either:
  - a character vector applied uniformly to all task outcomes, or
  - a named list with any of `success`, `warning`, and `error` to enable conditional retention

---

## 4. Formal Contracts (Unambiguous)

### 4.1 Data Generator Contract

#### Signature
```r
data_generator <- function(data_spec, seed, task_ctx) {}
```

#### Design note
- Keep `data_generator` as a plain function in v1 for lightweight extension and serialization simplicity
- If generator metadata/capability introspection becomes necessary, add an optional wrapper later rather than promoting generators to S7 prematurely

#### Required return structure: `data_bundle`
```r
list(
  train = object,                      # required; structure is fitter-defined
  test = NULL | object,                # optional; structure is fitter-defined
  response = NULL | character,         # optional hint for standard supervised workflows
  true_params = named numeric,         # names exactly equal vars_of_interest
  vars_of_interest = character,        # non-empty unique
  references = NULL | named numeric,   # if supplied, names exactly equal vars_of_interest
  meta = named list()                  # optional scalar metadata only
)
```

#### Validation invariants
- `train` exists and is non-NULL
- `setequal(names(true_params), vars_of_interest)` is `TRUE`
- if `references != NULL`, `setequal(names(references), vars_of_interest)` is `TRUE`
- if `response != NULL`, it is a non-empty character vector or scalar label understood by the fitter
- No duplicate names in any named vector/list
- Structural validation of `train`/`test` is delegated to the fitter, because not all Bayesian workflows are univariate supervised `data.frame` problems

---

### 4.2 Fitter Contract (S7)

```r
Fitter <- S7::new_class(
  "Fitter",
  abstract = TRUE,
  properties = list(
    name = S7::new_property(S7::class_character),
    supports_predictions = S7::new_property(S7::class_logical, default = TRUE),
    supports_log_lik = S7::new_property(S7::class_logical, default = TRUE),
    supports_loo = S7::new_property(S7::class_logical, default = TRUE)
  ),
  methods = list(
    fit = function(data_bundle, fit_spec, seed, task_ctx) 
      S7::stop_method_not_implemented(),
    extract_draws = function(fit_result, variables = NULL, draw_ids = NULL, ...) 
      S7::stop_method_not_implemented(),
    predict = function(fit_result, newdata = NULL, seed = NULL, ndraws = NULL, draw_ids = NULL, ...) 
      S7::stop_method_not_implemented(),
    log_lik = function(fit_result, newdata = NULL) 
      S7::stop_method_not_implemented(),
    loo = function(fit_result) 
      S7::stop_method_not_implemented(),
    diagnostics = function(fit_result) 
      S7::stop_method_not_implemented()
  )
)
```

#### `fit()` return: `bayesim_fit_result` (S3 list)
```r
structure(
  list(
    success = TRUE/FALSE,              # scalar logical
    fit = NULL | any,                  # backend object; may be removed by retention
    draws = NULL | matrix,             # S x P, colnames required
    diagnostics = named list,          # scalar values or named numeric vectors
    timing = list(total = numeric(1), warmup = numeric(1), sample = numeric(1)),
    warnings = character(),            # warning messages captured during fit
    error = NULL | condition
  ),
  class = "bayesim_fit_result"
)
```

#### Invariants
- `success == FALSE` requires non-NULL `error`
- `success == TRUE` requires `error == NULL`
- `timing$total >= 0`
- If `draws != NULL`, matrix must have column names

---

### 4.3 Metric Contract (S7)

```r
Metric <- S7::new_class(
  "Metric",
  properties = list(
    name = S7::new_property(S7::class_character),
    needs = S7::new_property(S7::class_character, default = character()),  # e.g. "predictions", "log_lik", "loo"
    required = S7::new_property(S7::class_logical, default = FALSE)        # if TRUE and fails -> task failure
  ),
  methods = list(
    compute = function(fit_result, data_bundle, context, task_ctx) 
      S7::stop_method_not_implemented()
  )
)
```

#### Metric output schema
- Named list, non-empty names
- Values must be one of:
  - scalar atomic (`logical`, `integer`, `double`, `character`)
  - short fixed-width named numeric vector
- No nested data frames/matrices in final metric output row
- Engine flattens only fixed-width outputs with prefix `<metric_name>__<field>`
- Per-observation, per-draw, or high-cardinality per-parameter outputs must be written to an artifact file or a long-format side table, with the result row storing only a pointer/hash plus compact summaries

---

## 5. Execution Model and Determinism

### 5.1 Task identity and ordering
Each task has deterministic key:
```r
task_id = sprintf(
  paste0("d%0", data_width, "d_f%0", fit_width, "d_r%0", rep_width, "d"),
  data_idx,
  fit_idx,
  rep_idx
)
```
The global task table is created once, sorted by integer task columns (`data_idx`, `fit_idx`, `rep_idx`), then persisted with a canonical string `task_id` derived from those columns.

If `task_grid` is supplied by the user, canonicalize and persist it directly; otherwise build the task table from the Cartesian product of `data_grid`, `fit_grid`, and `n_replicates`.

### 5.2 RNG policy (mandatory)
- Set RNG kind once at run start:
```r
RNGkind("L'Ecuyer-CMRG")
set.seed(config$seed)
```
- Precompute one independent RNG stream per task using `parallel::nextRNGStream()` (not by consuming random draws)
- Store `rng_seed` in task table
- In worker, before any random call:
```r
.Random.seed <- task$rng_seed
```
- Pass explicit seeds through to backends that maintain their own RNG state (for example Stan/brms backends)
- Use `future_lapply(..., future.seed = FALSE)` to avoid reseeding by backend
- **Never derive randomness from wall time, PID, or worker index**

### 5.3 Determinism guarantees
Same `config_fingerprint` + same package versions + same backend + same task table + same thread settings => identical task assignment, seed assignment, and sorted result ordering.

For numerically deterministic backends this should yield identical outputs; for backends with platform/thread/library variability, the guarantee is reproducibility within documented backend limits and tolerance-based equivalence where necessary.

---

## 6. Error Handling Strategy (Complete)

### 6.1 Condition classes
Use `rlang::abort()` with classes:

**Fatal (stop run):**
- `bayesim_config_error`
- `bayesim_contract_error`
- `bayesim_checkpoint_error`
- `bayesim_internal_error`

**Task-level recoverable:**
- `bayesim_data_error`
- `bayesim_fit_error`
- `bayesim_metric_error`

### 6.2 Policy
- Data or fit failure marks task `failed`; simulation continues unless `failed_count > max_errors`
- Metric failure:
  - If metric `required = FALSE`: emit NA fields + `<metric>__error_class`, continue
  - If metric `required = TRUE`: task becomes `failed`
- All errors captured with:
  - class, message, call (optional), and trimmed traceback string
- Ledger must record one terminal status per task: `success|failed|skipped`

---

## 7. Checkpoint/Resume Protocol (Final)

### 7.1 Directory layout
```
result_path/
├── run_manifest.json           # run-level metadata and schema versions
├── latest.json                 # pointer to latest valid checkpoint ID
└── checkpoints/
    └── cp_000001/
        ├── meta.json
        ├── ledger.rds          # authoritative checkpoint ledger
        ├── results.rds         # authoritative terminal task rows
        ├── result_chunks/      # optional columnar chunks for large runs
        ├── artifacts/          # optional task artifacts when retained
        └── checksums.json
```

### 7.1a Storage backend policy
- The authoritative v1 checkpoint/resume path is `RDS` + JSON because it is CRAN-safe, simple, and does not make core functionality depend on a suggested package
- Optional columnar persistence/export may be added behind `checkpoint_format = "parquet"` after benchmarking
- If a columnar backend is enabled, record the backend name and version in the manifest
- Benchmark candidates may include `r-polars` and `arrow`, but neither should become a correctness dependency without clear wins in write throughput, read performance, file size, and installation/CRAN stability

### 7.2 Schema versions
- `run_schema_version`: increments on any incompatible on-disk format change
- `result_schema_version`: increments on result-column contract changes
- Resume only allowed when both versions are compatible

### 7.3 Atomic write protocol
1. Write checkpoint to `cp_XXXXXX.tmp/` in the same parent directory as the target checkpoint
2. Flush files and write checksums
3. Validate tmp checkpoint (read-back + row-count + checksum)
4. Rename tmp dir to final `cp_XXXXXX/` (single filesystem operation on the same filesystem)
5. Atomically replace `latest.json` by write-temp-then-rename

If rename fails because the destination is on a different device or filesystem, abort with a checkpoint error by default; an explicit non-atomic fallback may be added later, but it must be opt-in and documented as weaker than the local-filesystem guarantee.

### 7.4 Config fingerprint
Fingerprint is hash of normalized `config_spec` excluding runtime-only fields:

**Excluded:** `result_path`, `checkpoint_every`, `progress`

**Included:** data/fit grids, generator/fitter identifiers + versions, metrics, retention, seed, n_replicates, checkpoint format, backend options, backend thread settings

For custom functions/closures, fingerprint a normalized representation of the function body plus identifiable package/namespace provenance. If a function cannot be normalized safely (for example because it captures opaque environment state), resume must be rejected and the user must start a fresh run with `run_simulation(..., resume = "never")`.

Manifest-only resume is allowed only when the fitter/data generator/metrics are rehydratable from the stored identifiers and versions; otherwise the user must supply `config` to `resume_simulation()`.

Resume requires an exact compatible match; otherwise the user must start a fresh run with `run_simulation(..., resume = "never")`.

### 7.5 Resume algorithm
1. Load `latest.json`; scan backward for most recent valid checkpoint
2. Validate checkpoint integrity and fingerprint
3. Rehydrate executable components from the manifest when possible; otherwise require user-supplied `config`
4. Rebuild full task table deterministically
5. Mark tasks terminal from checkpoint ledger
6. Execute only pending tasks
7. Merge prior + new results by `task_id`; treat unexpected duplicate terminal rows as an integrity error unless they are byte-identical or explicitly marked as a safe recovery rewrite
8. Write new checkpoint and final result object

---

## 8. Memory Management Strategy

### 8.1 Retention levels
- `minimal`: metrics only
- `standard` (default): metrics + diagnostics
- `debug`: metrics + diagnostics + draws + fit + data

User-facing API may accept either named profile or explicit set from:
```r
c("metrics", "diagnostics", "draws", "predictions", "fit", "data", "warnings")
```

### 8.2 Runtime controls
- Process tasks in chunks (`chunk_size` internal; default 50)
- After each chunk:
  - Append checkpoint
  - Drop large intermediates
- Call `gc()` only on the controller process after chunk commit or under clear memory pressure; never inside the tight per-task worker loop
- **Never keep all task-level fit objects in memory**
- If `retain` includes large artifacts, write per-task artifact files and keep only path in result row

### 8.3 Result row size guardrail
If serialized row size exceeds threshold (e.g. 5MB):
- Move payload to artifact file
- Replace in-row value with pointer + hash
- Log warning once per run

---

## 9. Testing Strategy (Including Stress)

### 9.1 Unit tests (contracts)
- Validators for `data_bundle`, `fit_result`, metric outputs
- S7 interface conformance tests for custom fitter stubs
- Error class and policy tests (fatal vs recoverable)

### 9.2 Determinism tests
- Same seed, sequential vs multisession => identical sorted results
- Same seed, interrupted+resume vs single uninterrupted run => identical
- Different seeds => at least one metric differs
- Backend-seed propagation tests for Stan/brms adapters

### 9.3 Checkpoint tests
- Atomicity under simulated crash during write (tmp dir remains, latest unchanged)
- Corruption recovery (invalid latest checkpoint falls back to previous valid)
- Fingerprint mismatch rejects resume and instructs the user to start a fresh run with `resume = "never"`
- Duplicate terminal row detection fails loudly
- Cross-filesystem checkpoint path is rejected with a clear error unless a weaker fallback mode is explicitly enabled

### 9.4 Stress/performance tests (skip on CRAN)
- 10,000 lightweight tasks with mock fitter under `retain = "minimal"`
- Memory growth budget: bounded and sublinear vs task count
- Long-run interruption test with multiple resume cycles
- High error-rate scenario to confirm graceful degradation and `max_errors` handling

### 9.5 Compatibility tests
- brms fitter with `cmdstanr` as the default backend on smoke-size workloads
- brms fitter with `rstan` as an optional/configurable backend where available
- Custom fitter fixture not importing brms to verify backend-agnostic design

### 9.6 Storage backend benchmarks
- Benchmark authoritative `RDS` checkpoints against optional columnar sinks on representative workloads
- Compare write time, resume latency, file size, dependency burden, and CRAN/platform reliability
- Do not promote `parquet` to the default resume path unless it wins clearly enough to justify the extra dependency surface

---

## 10. Implementation Roadmap (Realistic)

### Phase 0: Lock Persistence + Reproducibility Contract
- Rewrite the reproducibility guarantee to be backend-qualified rather than universally bitwise
- Lock the authoritative checkpoint format (`RDS` + JSON) and define the optional columnar benchmark track
- Define fingerprint policy for custom closures and backend thread settings
- **Deliverable:** Persistence/reproducibility decisions are explicit and testable

### Phase 1: Contracts and Skeleton (Day 1-2)
- Implement S7 classes (`SimulationConfig`, `Fitter`, `Metric`)
- Implement validators and S3 result constructors
- Add error classes and normalization utilities
- **Deliverable:** Contract tests all passing

### Phase 2: Deterministic Engine Core (Day 3-4)
- Build task table + RNG stream assignment with `parallel::nextRNGStream()`
- Implement worker execution and metric evaluation with per-task error capture
- **Deliverable:** Deterministic unit/integration tests passing sequentially

### Phase 3: Checkpoint + Resume (Day 5-6)
- Implement checkpoint writer/reader, checksum validation, latest pointer
- Implement resume merge logic and fingerprint checks
- **Deliverable:** Interruption/resume equivalence tests passing

### Phase 4: Memory/Retention and Artifact Externalization (Day 7)
- Implement retention profiles and row-size guardrail
- Add chunked processing and controller-side GC strategy
- **Deliverable:** Stress test with bounded memory in CI-optional profile

### Phase 5: brms Adapter and Docs (Day 8-9)
- Implement `BrmsFitter` using final interface with `cmdstanr` as default backend and `rstan` as configurable alternative
- Add vignettes for custom fitter/metric and reproducibility guarantees
- **Deliverable:** End-to-end smoke run + documented extension examples

### Phase 6: Optional Columnar Backend Benchmark
- Benchmark `r-polars` and `arrow` against the baseline `RDS` checkpoint path
- Promote a columnar backend only if it materially improves throughput/size without compromising portability or CRAN usability
- **Deliverable:** Storage recommendation backed by numbers rather than preference

---

## 11. Code Patterns (Reference Snippets)

### Safe task wrapper
```r
run_task_safe <- function(task, config_spec, fitter, metrics) {
  tryCatch(
    run_task_impl(task, config_spec, fitter, metrics),
    error = function(e) {
      list(
        task_id = task$task_id,
        status = "failed",
        error_class = class(e)[1],
        error_message = conditionMessage(e),
        metrics = NULL
      )
    }
  )
}
```

### Metric failure downgrade (non-required metric)
```r
compute_metric_safe <- function(metric, fit_result, data_bundle, context, task_ctx) {
  tryCatch(
    metric@compute(fit_result, data_bundle, context, task_ctx),
    error = function(e) {
      if (isTRUE(metric@required)) {
        rlang::abort(conditionMessage(e), class = "bayesim_metric_error")
      }
      list(
        value = NA_real_,
        error_class = class(e)[1],
        error_message = conditionMessage(e)
      )
    }
  )
}
```

### Atomic JSON write
```r
write_json_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE)
  ok <- file.rename(tmp, path)
  if (!ok) {
    rlang::abort("Atomic rename failed", class = "bayesim_checkpoint_error")
  }
}
```

### RNG stream precomputation
```r
create_task_rng_streams <- function(global_seed, n_tasks) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(global_seed)
  
  streams <- vector("list", n_tasks)
  seed <- .Random.seed
  for (i in seq_len(n_tasks)) {
    streams[[i]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }
  
  streams
}
```

### Backend seed propagation
```r
prepare_backend_seed <- function(task_seed, backend = c("generic", "cmdstanr", "rstan")) {
  backend <- match.arg(backend)
  .Random.seed <- task_seed

  if (backend == "generic") {
    return(list(r_seed = task_seed))
  }

  list(
    r_seed = task_seed,
    backend_seed = as.integer(abs(task_seed[2]) %% .Machine$integer.max)
  )
}
```

---

## 12. Definition of Done

- [ ] `devtools::check()` passes with no errors/warnings
- [ ] Contract, determinism, checkpoint, and stress tests pass in CI
- [ ] Resume behavior verified after forced interruption
- [ ] Custom fitter example works without brms dependency
- [ ] Reproducibility guarantees are documented with backend/platform caveats instead of universal bitwise claims
- [ ] Storage backend decision is benchmarked and justified
- [ ] Documentation includes:
  - [ ] Contract reference
  - [ ] Reproducibility policy
  - [ ] Checkpoint/resume behavior
  - [ ] Retention/memory tuning guide

---

## 13. Risks and Watch-Outs

| Risk | Mitigation |
|------|------------|
| Hidden nondeterminism from backend internals | Fixed RNG streams, explicit backend seeds, documented backend-qualified guarantees |
| Checkpoint bloat from retained artifacts | Retention defaults and row-size externalization |
| Wrong storage format choice for core persistence | Keep `RDS` authoritative first; benchmark columnar backends before promoting them |
| Overly complex abstractions | Keep exactly one extension path (`Fitter` + `Metric`), avoid extra plugin layers |
| Parallel serialization of fitter objects | Precompiled models are NOT serialized; workers compile locally or receive minimal specs |
| Stan backend differences | Default brms to `cmdstanr`, support `rstan` as configurable, and test both where available |
| Resume false-positives with custom closures | Normalize and fingerprint function bodies/env provenance; reject ambiguous resumes |
| Cross-filesystem atomicity assumptions | Write temp dirs in-place and fail loudly when atomic rename cannot be guaranteed |

---

## 14. Escalation Triggers

**Only add complexity if:**
- Workloads exceed local filesystem limits (checkpoint size/time dominates runtime)
- Then introduce optional partitioned result sinks (for example Parquet via `r-polars`/`arrow`, Arrow dataset, or DB)
- **Do not add this until stress tests show current protocol is insufficient**

---

## Appendix A: Package Structure

```
bayesim/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── bayesim-package.R
│   ├── errors.R                    # Error condition classes
│   ├── contracts.R                 # S7 interfaces + validators
│   ├── simulation-config.R         # S7 SimulationConfig
│   ├── fitter.R                    # S7 Fitter abstract class
│   ├── metric.R                    # S7 Metric abstract class
│   ├── results.R                   # S3 result constructors
│   ├── task-grid.R                 # Task table generation
│   ├── rng.R                       # RNG stream management
│   ├── simulate.R                  # Main run_simulation()
│   ├── worker.R                    # Task execution
│   ├── checkpoint.R                # Checkpoint read/write
│   ├── resume.R                    # Resume logic
│   ├── retention.R                 # Memory management
│   ├── metric-registry.R           # Built-in metrics registry
│   ├── metrics-posterior.R
│   ├── metrics-diagnostics.R
│   ├── metrics-predictive.R
│   ├── brms-fitter.R               # BrmsFitter implementation
│   └── utils.R
│
├── tests/
│   ├── testthat.R
│   └── testthat/
│       ├── helper-fixtures.R
│       ├── helper-mocks.R
│       ├── test-contracts.R
│       ├── test-validators.R
│       ├── test-errors.R
│       ├── test-task-grid.R
│       ├── test-rng.R
│       ├── test-simulate.R
│       ├── test-checkpoint.R
│       ├── test-resume.R
│       ├── test-determinism.R
│       ├── test-retention.R
│       ├── test-metrics.R
│       ├── test-brms-fitter.R
│       └── test-stress.R
│
├── vignettes/
│   ├── getting-started.Rmd
│   ├── custom-fitters.Rmd
│   ├── custom-metrics.Rmd
│   ├── reproducibility.Rmd
│   └── memory-management.Rmd
│
└── man/
```

---

## Appendix B: DESCRIPTION

```yaml
Package: bayesim
Title: Simulation Framework for Bayesian Modeling
Version: 1.0.0

Authors@R: c(
    person("Maximilian", "Scholz", , "research.scholz@mailbox.org", 
           role = c("aut", "cre")),
    person("Paul-Christian", "Bürkner", , "paul.buerkner@gmail.com", 
           role = "aut")
  )

Description: A modern simulation framework for Bayesian modeling studies.
    Provides extensible tools for running complex simulation studies with
    deterministic reproducibility, checkpoint/resume capabilities, and 
    memory-bounded execution. Users can extend with custom fitters, metrics,
    and data generators.

License: GPL (>= 3)

Config/testthat/edition: 3
Encoding: UTF-8
Roxygen: list(markdown = TRUE)

Imports:
    S7,
    cli (>= 3.6.0),
    digest (>= 0.6.30),
    dplyr (>= 1.1.0),
    future.apply (>= 1.11.0),
    jsonlite (>= 1.8.0),
    lifecycle (>= 1.0.3),
    matrixStats (>= 1.0.0),
    posterior (>= 1.4.0),
    rlang (>= 1.1.0),
    stats,
    tibble (>= 3.2.0),
    withr (>= 2.5.0)

Suggests:
    arrow (>= 14.0.0),
    brms (>= 2.22.0),
    cmdstanr (>= 0.7.0),
    rstan (>= 2.32.0),
    testthat (>= 3.1.0),
    knitr,
    rmarkdown

VignetteBuilder: knitr
```
