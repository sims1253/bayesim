# bayesim 2.0 — Remaining Review Backlog

Updated: 2026-07-11, branch `feat/plan-implementation` after commits
`0c2f176` and `cf3ddf3`.

This document contains only work that remains useful after the July 2026
review-hardening pass. Resolved findings were removed rather than retained as a
second changelog; implementation history belongs in `NEWS.md` and git.

## Current assessment

The confirmed correctness, retention, RNG, model-bank, documentation, and
repository-hygiene findings from the original review are fixed. The default
test suite, formatting, lint, and package check pass. The remaining work is
mostly scalability, API depth, CI cost, and a few robustness edges rather than
evidence that the core simulation results are wrong.

## 1. Correctness and robustness edges

### 1.1 Make summary grouping collision-safe

`summarize_simulation()` currently constructs group identifiers by pasting
condition values together with `"|"`. Distinct cells can collide when a value
contains the delimiter (for example `c("a|b", "c")` versus
`c("a", "b|c")`), and a literal `"NA"` may be indistinguishable from a missing
value. Replace the pasted key with collision-safe row grouping, such as
an explicitly NA-safe stable row-identity helper, and add tests for delimiters,
missing values, factors, and mixed condition types.

### 1.2 Clarify adaptive-stop progress reporting

The stopping rule now correctly requires the MCSE target in every condition
cell. Its console message still says the threshold was reached after
`n_completed` "reps", although that value is the total number of completed
tasks across all cells. For multi-condition studies, report "completed tasks"
or show the minimum/maximum replicate count per condition.

### 1.3 Strengthen precompiled-brms prior detection

The model bank warns when a model row supplies no explicit prior. A partially
specified prior can still leave data-dependent brms defaults embedded from the
template dataset. If practical, compare the effective/default prior set and
warn whenever a data-dependent prior class remains implicit. Until then, keep
the documented rule stronger than the detection: precompiled studies should
specify every relevant prior explicitly.

## 2. Checkpoint scalability

The number of checkpoint snapshots is now bounded by `keep_checkpoints`
(default 2), and resumed results are cached so each batch does not re-read and
re-hash the prior checkpoint. However, every checkpoint is still a complete
cumulative snapshot. The engine therefore rewrites the full task ledger and an
ever-growing results table after each batch; with a fixed batch size, total
bytes written grow quadratically with the number of batches.

Before redesigning the format, benchmark realistic 10k, 100k, and 1M-task
studies and record controller time, bytes written, resume time, and peak memory.
If write amplification is material, move to append-only per-batch result files
plus periodic compaction. Preserve the properties users already rely on:

- atomic visibility of a completed checkpoint;
- checksums and corruption fallback;
- deterministic resume and task de-duplication;
- RDS checkpoint behavior and the optional parquet summary sidecar;
- bounded retention of historical snapshots.

Treat the checkpoint settings that now travel through configuration,
execution, and storage as a candidate `checkpoint_policy` value object during
that redesign, rather than adding more primitive parameters to each layer.

## 3. Performance follow-up

The execution path now resolves task rows once and updates batch statuses
vectorially. Validate the improvement with a controller-only benchmark at
100k–1M tasks; tests at ordinary package scale will not expose remaining
allocation or serialization costs.

Several older internal helpers (`setup_global_rng()`, ID-based
`get_task_spec()`, and `update_task_status()`) no longer serve the execution
path and are exercised mainly by their own tests. Audit whether they are useful
extension seams. If not, remove them and their generated documentation before
they become accidental compatibility commitments.

## 4. Statistical extensions

### 4.1 Implement a calibrated multi-chain SBC band

For `L > 1`, `adjust_gamma()` explicitly falls back to the exact `L = 1` band
and informs the user that it is conservative. Port or independently implement
the correlated/multi-chain calibration from the reference method if multi-chain
rank ECDFs become a supported first-class workflow. Add numerical parity tests
against a trusted implementation, not only shape/property tests.

### 4.2 Decide whether discrete-parameter SBC is in scope

`rank_metric()` documents that its strict-below rank is intended for continuous
posteriors and recommends randomized tie-breaking for discrete or boundary-mass
parameters. If bayesim expands into those models, add an explicit tie strategy
(including deterministic seeding and tests). Otherwise keep the limitation
prominent and avoid implying universal SBC support.

### 4.3 Broaden numerical regression tests selectively

The simultaneous-band implementation now has a fixed numerical regression
test. Additional high-value tests would compare several `(N, K, alpha)` cases
against the reference implementation and empirically verify global coverage.
Keep stochastic acceptance tolerances predeclared so the tests do not become
flaky.

## 5. API and module design

### 5.1 Add a metric-column accessor

Flattened names such as `posterior_summary__q_lower__x` are stable but awkward
for interactive analysis. A small exported helper such as
`metric_cols(result, "posterior_summary")`, optionally returning a named
selection or tidy long form, would reduce stringly typed user code. Design it
around the existing flattening contract rather than adding a second naming
scheme.

### 5.2 Separate study identity from runtime policy in a future breaking release

`simulation_config()` remains broad because it combines the statistical study
definition with retention, checkpoint, error, daemon, and stopping policy. The
fingerprint already distinguishes those concepts. In the next deliberate API
break, consider a `runtime_options()` object or arguments on `run_simulation()`.
Do not change this piecemeal: checkpoint/resume compatibility, printing,
manifests, targets integration, and documentation should move together.

### 5.3 Deepen retention and analysis seams when they next change

`build_metric_context()` currently understands both character and
context-specific list retention forms. Move that knowledge behind a retention
predicate (for example `retention_may_request(retain, "predictions")`) when the
retention API next evolves.

`R/analysis.R` contains generic aggregation, SBC extraction, and plotting. It
also performs separate per-group traversals when building SBC ECDF and band
data. If this layer grows, extract a per-group SBC helper and split aggregation,
SBC accessors, and plotting into focused modules. This is maintainability work,
not a reason for immediate churn.

## 6. CI and test feedback

The expensive Stan-backed SBC acceptance test is opt-in locally but enabled on
the full five-platform R CMD check matrix. That gives strong coverage at a high
feedback and infrastructure cost. Add a fast no-Stan test job as the first PR
gate, then consider running the acceptance suite on one Linux job (plus a
scheduled or main-branch run) while retaining cross-platform package checks.

The mirai parity test can warn that `'package:bayesim' may not be available
when loading` under `devtools::load_all()`. Reproduce the test against an
installed package. If the warning disappears, adjust the development test
harness or suppress only that verified load-all artifact; if it persists,
fix worker serialization/package loading before treating it as harmless.

Add lightweight execution checks for examples currently marked `eval = FALSE`
where feasible, especially the targets pipeline. External-package API drift is
otherwise discovered only by manual vignette review.

## 7. Scope and release readiness

The package remains lifecycle-experimental. Before declaring a stable release,
set explicit maturity criteria for the secondary surfaces (`stop_on`, parquet
sidecars, `report()`, and `rstar_metric()`): supported inputs, performance
expectations, failure semantics, and minimum tests. Prefer graduating or
removing weak surfaces over adding more adjacent features.

Before any CRAN submission, perform a dedicated policy/dependency audit rather
than assuming the current GitHub-oriented setup is ready. In particular, decide
whether the remaining `Remotes: stan-dev/cmdstanr` entry is still necessary,
verify examples/vignettes without optional toolchains, and make the intended
release version agree across `DESCRIPTION`, `NEWS.md`, and the tag.

## Suggested order

1. Fix collision-safe grouping and adaptive-stop messaging.
2. Benchmark large task grids and cumulative checkpoint writes.
3. Reproduce the mirai warning with an installed package and add a fast CI gate.
4. Decide checkpoint-format and runtime-policy redesign together.
5. Add API/statistical extensions only when concrete user studies require them.
