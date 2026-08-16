# bayesim Improvement Plan

## Objective

Strengthen bayesim’s run lifecycle, statistical correctness, extension
interfaces, usability, and test strategy without expanding the suite
into a large collection of low-value unit and smoke tests.

The target is a package whose main scientific workflows are verified
through a small number of end-to-end integration tests, whose checkpoint
and resume behavior is safe and predictable, and whose public extension
interfaces do not require knowledge of internal result structures.

## Guiding principles

- Treat
  [`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
  and the objects returned through its public interface as the primary
  test surface.
- Prefer deep modules: small interfaces hiding substantial behavior.
- Separate scientific study identity from runtime and storage policy.
- Never silently overwrite or mix results from incompatible studies.
- Persist canonical task outcomes that round-trip without reconstructive
  guessing.
- Distinguish fixed-truth performance studies from varying-truth
  calibration studies.
- Replace superseded internal tests when integration coverage exists; do
  not merely add another layer.
- Keep Stan-backed validation out of the fast pull-request feedback
  loop.

## Priority 0: Correctness and data safety

### 1. Define lifecycle invariants

Before restructuring code, document and agree on:

- What happens when `result_path` does not exist, is empty, contains a
  matching run, contains an incompatible run, or contains corrupt state.
- The exact meanings of `resume = "auto"`, `"never"`, and `"must"`.
- Which task states are terminal and which remain eligible for resume.
- Whether failed tasks are retried, and under what explicit policy.
- What changing retention or adaptive-stop settings on resume means.
- Which configurations support `resume_simulation(path)` without
  supplying the original configuration.
- Whether truth must be fixed within a performance-measure condition
  cell.

Recommended path behavior:

- A matching fingerprint may resume.
- An incompatible or ambiguous directory aborts with a diagnostic.
- `resume = "never"` requires an empty/new directory unless an explicit
  overwrite option is supplied.
- Checkpoints from different fingerprints are never stored in the same
  run directory.

### 2. Add black-box regression tests for known defects

Add public-interface tests for:

- Incompatible configuration reuse of an existing `result_path`.
- `max_errors = 0` and other error-budget edge cases.
- Continuation after a policy-stopped run.
- Resume parity for summary rows and task-result structure.
- Retention behavior across resume.
- Varying-truth bias MCSE.
- Invalid interval probabilities.
- Duplicate custom task identities.
- A sequential/parallel prediction workflow with non-`NA` metric values.

These tests should be written before changing behavior and should assert
observable outcomes rather than internal helper calls.

### 3. Fix the task state machine

Replace the current overloaded use of `"skipped"` with explicit states
or stop metadata. A suitable state model is:

- `pending`
- `running`
- `success`
- `failed`
- `not_run_policy_stop`
- `cancelled`

Only genuinely completed tasks should be terminal during resume.
Unexecuted tasks left by `max_errors` or adaptive stopping must remain
eligible for later execution.

Completion criteria:

- Raising the error budget and resuming executes previously unexecuted
  tasks.
- A stricter adaptive target can resume a previously stopped run.
- `max_errors = 0` has documented and tested semantics.
- Interrupted and uninterrupted runs produce equivalent canonical
  outcomes.

### 4. Correct performance-measure semantics

Decide between two supported modes:

1.  Require truth to be constant within each condition cell for
    Morris-style bias, empirical SE, and model-SE comparisons; or
2.  Add separately named varying-truth measures with appropriate
    formulas and interpretations.

At minimum:

- Bias MCSE must be based on the replicate-level error when truth
  varies.
- Empirical SE must not be presented with a fixed-truth interpretation
  when truth varies.
- SBC documentation must not imply that all Morris measures remain
  directly interpretable under prior-predictive varying truth.
- Statistical formulas need fixed reference-fixture tests.

## Priority 1: Architectural restructuring

### 5. Separate study identity from run policy

Introduce two internal modules while retaining
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
as a backward-compatible facade:

#### `StudySpec`

Contains only scientific identity:

- Data and fit grids, or an explicit task design.
- Data generator.
- Fitter.
- Metrics.
- Replicate count.
- Study seed.

It should be immutable after validation and should be the source of the
study fingerprint.

#### `RunPolicy`

Contains operational choices:

- Worker configuration.
- Result path and storage format.
- Checkpoint cadence and retention.
- Error budget and retry policy.
- Adaptive stopping.
- Progress and verbosity.

Completion criteria:

- Fingerprinting consumes only `StudySpec`.
- Resume compatibility decisions explicitly account for policies that
  affect persisted artifacts or task completion.
- Core run code no longer mutates the scientific configuration.

### 6. Introduce a run-store seam

Create an internal `RunStore` interface with at least two real adapters:

- In-memory adapter for runs without `result_path` and for tests.
- Filesystem adapter for durable checkpointed runs.

The store should own:

- Run-directory initialization and collision checks.
- Manifest and fingerprint validation.
- Task ledger persistence.
- Appending canonical task outcomes.
- Reading outcomes for resume and final result assembly.
- Integrity verification and corruption fallback.

Avoid rewriting the complete accumulated result set after every batch.
Prefer append-only result shards plus an atomically updated ledger, with
optional periodic compaction.

Completion criteria:

- Checkpoint work grows approximately linearly with completed tasks.
- Lightweight run memory does not grow through retained copies of every
  full task result.
- The in-memory and filesystem adapters pass the same lifecycle
  integration tests.
- Corrupt latest state falls back safely without overwriting prior
  history.

### 7. Persist a canonical `TaskOutcome`

Replace lossy reconstruction from flattened summary columns with a
canonical persisted outcome containing distinct fields for:

- Task identity and status.
- Metrics.
- Diagnostics.
- Truth.
- Warnings.
- Error information.
- Timing.
- Retained or externalized artifacts.

The flat summary should be a derived view, not the canonical persisted
object.

Completion criteria:

- A task outcome is structurally equivalent before and after checkpoint
  round-trip.
- Resume does not move truth or diagnostics into `$metrics`.
- Artifact availability is consistent and explicitly reported.
- Retention changes on resume either have well-defined limitations or
  are rejected when they cannot be fulfilled.

### 8. Redesign adaptive scheduling

Schedule adaptive studies by replicate rounds across all condition
cells, rather than completing one condition before moving to the next.

Replace exact modulo checks with threshold-crossing checks. Persist the
next scheduled check or the number of completed replicate rounds.

Completion criteria:

- Every condition cell reaches `min_reps` before any cell receives
  excessive additional work.
- Checks are not missed when checkpoint and adaptive intervals differ.
- Stop reasons and the precision state that triggered them are
  persisted.

## Priority 2: Extension interfaces and metric schema

### 9. Deepen the fitter interface

Reduce the knowledge required to implement a fitter:

- Export a supported fit-result constructor or replace the manually
  assembled S3 result with a public result type.
- Require only core fitting behavior.
- Supply default unsupported implementations for optional capabilities.
- Make capability declarations and method availability consistent.
- Provide a conformance function that accepts a user-supplied data
  bundle and fit specification instead of assuming `y ~ x`.
- Validate prediction, draw, log-likelihood, and diagnostic outputs at
  the task seam on representative/first use.

Contract violations such as prediction-length or matrix-orientation
mismatches should fail clearly; metrics should not truncate or recycle
invalid outputs.

Completion criteria:

- A package-external fitter can be implemented using only exported
  names.
- A fitter without prediction or LOO support does not need dummy
  methods.
- The custom-fitter vignette is executable and exercises every
  capability it claims to demonstrate.

### 10. Add field-level metric metadata

Replace one `summary_type` per metric with metadata for each emitted
field. Useful metadata includes:

- Field role: estimate, binary outcome, count, diagnostic, rank,
  artifact.
- Aggregation rule and MCSE method.
- Interval probability or nominal level.
- Parameter dimension.
- Units or scale where relevant.
- Externalization policy.

Validate metric names, `needs`, `required`, output schemas, and
constructor parameters such as probabilities at construction time.

Completion criteria:

- Binary per-parameter coverage and across-parameter means are not
  assigned the same Bernoulli MCSE blindly.
- Counts such as `n_obs` are not accidentally treated as scientific
  measures.
- [`plot_coverage()`](https://sims1253.github.io/bayesim/reference/plot_coverage.md)
  obtains or validates the correct nominal interval level.
- Results retain enough metadata for downstream analysis without name
  heuristics.

### 11. Improve analysis and plotting behavior

- Make grouped coverage plots visually distinguish or facet condition
  cells.
- Validate all requested grouping and plotting columns.
- Reject degenerate R-squared cases or return an explicit undefined
  result.
- Ensure failures, skipped tasks, and missing values have consistent
  inclusion rules in summaries.
- Keep a tidy long-form view available for metrics while retaining a
  convenient wide per-task summary.

## Priority 3: Usability and documentation

### 12. Establish three golden workflows

Make documentation center on three fully executable paths:

1.  Fixed-truth estimator-performance study.
2.  Prior-predictive SBC/calibration study.
3.  Package-external custom fitter and metric.

Every example must use only exported interfaces and must produce
meaningful, non-`NA` results for the behavior it demonstrates.

### 13. Clarify resume and reproducibility limitations

Document:

- When configless resume works.
- Why inline/global functions generally require supplying the
  configuration again.
- Which closure values participate in fingerprinting.
- What backend, package-version, and platform changes mean for
  reproducibility.
- What retention changes can and cannot recover from prior checkpoints.

### 14. Simplify common interaction

- Keep
  [`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
  compatible, but present scientific arguments separately from
  execution/storage options in documentation.
- Provide a clear run summary with paths, completed/pending counts,
  failures, and resume instructions.
- Make verbosity independent from progress display.
- Improve discovery of metric columns through documented selectors and
  tidy views rather than requiring users to memorize double-underscore
  names.

## Integration-focused testing strategy

### Fast suite: every pull request

Keep fewer than roughly 20 substantive tests across these workflows:

1.  **Complete analytic study**
    - Two data conditions and at least one held-out test set.
    - Posterior summaries, prediction metrics, coverage, performance
      measures, and analysis output.
2.  **Lifecycle parity**
    - Uninterrupted versus checkpoint/resume.
    - In-memory versus filesystem store.
    - Sequential versus `workers = 2`.
    - Compare canonical outcomes excluding wall-clock timing.
3.  **Failure and recovery**
    - Recoverable generator, fitter, and metric failures.
    - Fatal failure through parallel transport.
    - Error-budget stop followed by continuation.
    - Failure-reporting output.
4.  **External extension**
    - Custom generator, fitter, and metric defined outside bayesim.
    - Only exported interfaces allowed.
    - Prediction and log-likelihood shapes exercised, not merely
      constructed.
5.  **Adaptive multi-condition study**
    - Replicate-round scheduling.
    - Precision threshold reached across all cells.
    - Resume under a stricter threshold.
6.  **Analytic SBC acceptance**
    - Exact conjugate prior-predictive workflow.
    - Rank bounds, retained rank counts, deterministic summaries, and
      stable calibration statistics.

### Focused lower-level tests

Retain targeted tests for behavior that is difficult to diagnose through
the public workflow:

- Statistical formula parity against fixed fixtures.
- RNG isolation and stream derivation.
- Fingerprint canonicalization.
- Checksum and corruption fallback.
- Schema migration fixtures.
- Parameterized fitter/metric conformance checks.

Once integration coverage exists, remove repetitive tests for
constructor defaults, individual class hierarchy entries, shallow
wrappers, and every minor invalid-input permutation.

### Backend suite: on merge or release

Run on one cached Linux job:

- One tiny `BrmsFitter` workflow.
- One tiny `CmdStanFitter` workflow.
- Sequential/parallel parity.
- Draw, prediction, log-likelihood, LOO, and diagnostics conformance.
- An instrumented assertion of actual compilation count.
- Reference parity with `brms` and `loo` where scientifically important.

Do not install CmdStan on every operating-system/R-version package-check
job.

### CI layout

- **Pull request:** package check without CmdStan, fast integration
  suite, formatting, and linting.
- **Merge to master:** cached Stan backend suite and heavier statistical
  acceptance.
- **Release:** full supported R/OS matrix plus backend acceptance on the
  chosen reference platform.
- **Coverage:** fast/core suite only; use coverage to find untested
  behavior, not as a percentage target.
- Use explicit test tiers instead of `skip_on_cran()` as a speed-control
  mechanism.

## Suggested delivery sequence

### Milestone 1: Safety baseline

- Lifecycle invariants documented.
- Known-defect integration tests added.
- Path collision, `max_errors`, skipped/pending resume, and probability
  validation fixed.
- Varying-truth performance semantics corrected or rejected clearly.

### Milestone 2: Canonical lifecycle

- `StudySpec` and `RunPolicy` introduced internally.
- Task state machine implemented.
- Canonical `TaskOutcome` persisted and round-tripped.
- Resume parity achieved.

### Milestone 3: Scalable persistence and scheduling

- `RunStore` adapters implemented.
- Append-only/sharded result persistence replaces repeated full
  rewrites.
- Adaptive replicate-round scheduling implemented.
- Memory and checkpoint scaling benchmarks added.

### Milestone 4: Public extension redesign

- Supported fitter result constructor/interface exported.
- Optional fitter capabilities simplified.
- Runtime seam validation added.
- Field-level metric schemas introduced.

### Milestone 5: Usability and test replacement

- Three golden workflows documented and executable.
- Fast integration suite becomes the primary PR gate.
- Superseded internal tests removed.
- Backend checks moved to scheduled/release CI.
- Package-check, lint, format, and integration workflows all pass.

## Definition of done

The improvement program is complete when:

- Existing result data cannot be silently overwritten by an incompatible
  run.
- Every unexecuted task can be resumed under a compatible study
  definition.
- Checkpoint/resume is behaviorally equivalent to uninterrupted
  execution.
- Persisted task outcomes round-trip without field-role loss.
- Statistical summaries distinguish fixed- and varying-truth studies.
- Custom fitters and metrics use only documented exported interfaces.
- Fast integration tests exercise the complete analytic lifecycle with
  meaningful non-`NA` outputs.
- Stan integration is verified in a focused cached job rather than every
  fast package check.
- The test suite is smaller, more behavior-focused, and resilient to
  internal refactoring.
