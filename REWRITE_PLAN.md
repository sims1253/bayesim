# bayesim 1.0 Rewrite Plan (Final, Implementation-Ready)

**Version:** 2.0  
**Status:** Approved blueprint for implementation  
**Target:** CRAN-ready package with deterministic, resumable, memory-bounded simulation engine

---

## Bottom Line

This rewrite keeps the architecture simple: **S7 only for extension contracts** (`SimulationConfig`, `Fitter`, `Metric`), **S3 for all runtime and persisted data**.

The core engine must be deterministic across sequential/parallel/resume modes by assigning each task a precomputed L'Ecuyer RNG stream and never depending on scheduler order.

Checkpointing should be append-safe, atomic, and restartable using a validated ledger + chunked result files with strict schema versioning.

---

## Action Plan

1. Lock contracts first (interfaces, schemas, error classes, checkpoint format)
2. Build deterministic task execution + ledger before implementing advanced metrics
3. Add checkpoint/resume with corruption recovery and config fingerprint checks
4. Implement retention/memory controls and chunked writes
5. Finish with integration + stress tests (determinism, interruption, memory growth)

**Effort estimate:** Medium (1-2d) for contract scaffolding; Large (3d+) for full production implementation with robust testing.

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
run_simulation <- function(config, resume = FALSE, force_restart = FALSE, progress = TRUE) {}

# Config constructor
simulation_config <- function(
  data_grid,
  fit_grid,
  data_generator,      # function(data_spec, seed, task_ctx) -> data_bundle
  fitter,              # S7 Fitter instance
  metrics,             # character vector or list of Metric objects
  n_replicates = 1L,
  seed,
  result_path = NULL,
  checkpoint_every = 50L,
  retain = c("metrics", "diagnostics"),
  max_errors = Inf
) {}

# Registry
register_metric <- function(metric, overwrite = FALSE) {}
get_metric <- function(name) {}
list_metrics <- function() {}
unregister_metric <- function(name) {}
```

### Defaults
- `retain = c("metrics", "diagnostics")`
- `checkpoint_every = 50L`
- `max_errors = Inf` (user can fail-fast by setting finite threshold)

---

## 4. Formal Contracts (Unambiguous)

### 4.1 Data Generator Contract

#### Signature
```r
data_generator <- function(data_spec, seed, task_ctx) {}
```

#### Required return structure: `data_bundle`
```r
list(
  train = data.frame,                  # nrow >= 1
  test = NULL | data.frame,            # optional; if not NULL must contain required columns
  response = character(1),             # name of response column present in train/test
  true_params = named numeric,         # names exactly equal vars_of_interest
  vars_of_interest = character,        # non-empty unique
  references = named numeric,          # names exactly equal vars_of_interest
  meta = named list()                  # optional scalar metadata only
)
```

#### Validation invariants
- `setequal(names(true_params), vars_of_interest)` is `TRUE`
- `setequal(names(references), vars_of_interest)` is `TRUE`
- `response %in% names(train)` and if `test != NULL`, `response %in% names(test)`
- No duplicate names in any named vector/list
- No factor levels in `test` absent from `train` unless explicitly allowed by fitter

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
    extract_draws = function(fit_result, variables = NULL) 
      S7::stop_method_not_implemented(),
    predict = function(fit_result, newdata = NULL, seed = NULL) 
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
  - named numeric vector
- No nested data frames/matrices in final metric output row
- Engine flattens with prefix `<metric_name>__<field>`

---

## 5. Execution Model and Determinism

### 5.1 Task identity and ordering
Each task has deterministic key:
```r
task_id = sprintf("d%03d_f%03d_r%05d", data_idx, fit_idx, rep_idx)
```
The global task table is created once, sorted lexicographically by `task_id`, and persisted.

### 5.2 RNG policy (mandatory)
- Set RNG kind once at run start:
```r
RNGkind("L'Ecuyer-CMRG")
set.seed(config$seed)
```
- Precompute one RNG stream per task (integer vector seed state)
- Store `rng_seed` in task table
- In worker, before any random call:
```r
.Random.seed <- task$rng_seed
```
- Use `future_lapply(..., future.seed = FALSE)` to avoid reseeding by backend
- **Never derive randomness from wall time, PID, or worker index**

### 5.3 Determinism guarantees
Same `config_fingerprint` + same package versions + same task table => identical outputs (bitwise equality for numeric paths that are deterministic in backend).

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
        ├── ledger.parquet
        ├── results.parquet     # metrics + diagnostics, one row per terminal task
        ├── artifacts/          # optional task artifacts when retained
        └── checksums.json
```

### 7.2 Schema versions
- `run_schema_version`: increments on any incompatible on-disk format change
- `result_schema_version`: increments on result-column contract changes
- Resume only allowed when both versions are compatible

### 7.3 Atomic write protocol
1. Write checkpoint to `cp_XXXXXX.tmp/`
2. Flush files and write checksums
3. Validate tmp checkpoint (read-back + row-count + checksum)
4. Rename tmp dir to final `cp_XXXXXX/` (single filesystem operation)
5. Atomically replace `latest.json` by write-temp-then-rename

### 7.4 Config fingerprint
Fingerprint is hash of normalized `config_spec` excluding runtime-only fields:

**Excluded:** `result_path`, `checkpoint_every`, `progress`

**Included:** data/fit grids, generator/fitter identifiers + versions, metrics, retention, seed, n_replicates

Resume requires exact match unless `force_restart = TRUE`.

### 7.5 Resume algorithm
1. Load `latest.json`; scan backward for most recent valid checkpoint
2. Validate checkpoint integrity and fingerprint
3. Rebuild full task table deterministically
4. Mark tasks terminal from checkpoint ledger
5. Execute only pending tasks
6. Merge prior + new results by `task_id` (deduplicate by last-write-wins within same run_id)
7. Write new checkpoint and final result object

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
  - Call `gc()` when retained bytes exceed threshold or every N chunks
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

### 9.3 Checkpoint tests
- Atomicity under simulated crash during write (tmp dir remains, latest unchanged)
- Corruption recovery (invalid latest checkpoint falls back to previous valid)
- Fingerprint mismatch rejects resume unless `force_restart`

### 9.4 Stress/performance tests (skip on CRAN)
- 10,000 lightweight tasks with mock fitter under `retain = "minimal"`
- Memory growth budget: bounded and sublinear vs task count
- Long-run interruption test with multiple resume cycles
- High error-rate scenario to confirm graceful degradation and `max_errors` handling

### 9.5 Compatibility tests
- brms fitter (rstan/cmdstanr where available) on smoke-size workloads
- Custom fitter fixture not importing brms to verify backend-agnostic design

---

## 10. Implementation Roadmap (Realistic)

### Phase 1: Contracts and Skeleton (Day 1-2)
- Implement S7 classes (`SimulationConfig`, `Fitter`, `Metric`)
- Implement validators and S3 result constructors
- Add error classes and normalization utilities
- **Deliverable:** Contract tests all passing

### Phase 2: Deterministic Engine Core (Day 3-4)
- Build task table + RNG stream assignment
- Implement worker execution and metric evaluation with per-task error capture
- **Deliverable:** Deterministic unit/integration tests passing sequentially

### Phase 3: Checkpoint + Resume (Day 5-6)
- Implement checkpoint writer/reader, checksum validation, latest pointer
- Implement resume merge logic and fingerprint checks
- **Deliverable:** Interruption/resume equivalence tests passing

### Phase 4: Memory/Retention and Artifact Externalization (Day 7)
- Implement retention profiles and row-size guardrail
- Add chunked processing and periodic GC strategy
- **Deliverable:** Stress test with bounded memory in CI-optional profile

### Phase 5: brms Adapter and Docs (Day 8-9)
- Implement `BrmsFitter` using final interface
- Add vignettes for custom fitter/metric and reproducibility guarantees
- **Deliverable:** End-to-end smoke run + documented extension examples

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
  for (i in seq_len(n_tasks)) {
    streams[[i]] <- .Random.seed
    runif(1)  # advance stream
  }
  
  streams
}
```

---

## 12. Definition of Done

- [ ] `devtools::check()` passes with no errors/warnings
- [ ] Contract, determinism, checkpoint, and stress tests pass in CI
- [ ] Resume behavior verified after forced interruption
- [ ] Custom fitter example works without brms dependency
- [ ] Documentation includes:
  - [ ] Contract reference
  - [ ] Reproducibility policy
  - [ ] Checkpoint/resume behavior
  - [ ] Retention/memory tuning guide

---

## 13. Risks and Watch-Outs

| Risk | Mitigation |
|------|------------|
| Hidden nondeterminism from backend internals | Fixed RNG streams and deterministic task ordering |
| Checkpoint bloat from retained artifacts | Retention defaults and row-size externalization |
| Overly complex abstractions | Keep exactly one extension path (`Fitter` + `Metric`), avoid extra plugin layers |
| Parallel serialization of fitter objects | Precompiled models are NOT serialized; workers compile locally or receive minimal specs |
| Stan backend differences | Test both rstan and cmdstanr; document known differences |

---

## 14. Escalation Triggers

**Only add complexity if:**
- Workloads exceed local filesystem limits (checkpoint size/time dominates runtime)
- Then introduce optional partitioned result sinks (Arrow dataset or DB)
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
