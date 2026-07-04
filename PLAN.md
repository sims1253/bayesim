# bayesim 2.1 — Full-package review and redesign plan

**Status:** Approved work list from a whole-package ("thermonuclear") review, 2026-07-04.
**Supersedes:** the previous PLAN.md (F1–F8 remediation — implemented and verified on this branch).
**Audience:** an implementing agent with full access to this repository.
**Companion files:** `PLAN-bayesfam.md` and `PLAN-bayeshear.md` are standalone briefs for the
user's other repositories. Do NOT act on them here; they are hand-off documents.

---

## 0. Process requirements (read first)

- **Never run git commands that modify state** (add, commit, checkout, reset, stash, push...).
  The user handles all git operations. Read-only git (log, show, diff) is fine.
- **Never delete files without explicit user confirmation.** Flag deletion candidates instead.
- The package is pre-release (GitHub-only, lifecycle: experimental). **Breaking changes are
  allowed and expected** — do them cleanly (no deprecation shims for pre-release API), and
  record every one in NEWS.md under 2.0.0 (still unreleased; fold everything into that entry).
- After every workstream: `Rscript -e 'devtools::document()'` then
  `Rscript -e 'devtools::check(args = "--no-manual", error_on = "warning")'`.
  Green check is necessary but not sufficient — each item below names the behavior that must
  be tested. The SBC acceptance test (`tests/testthat/test-sbc-acceptance.R`) must stay green.
- Environment: R 4.6.x on Linux (WSL2), cmdstan via cmdstanr, mirai installed. Verify
  `packageVersion("purrr") >= "1.1.0"` and that `carrier` is available before starting C1.
- `README.Rmd` at the repo root was written by the reviewing agent **assuming this plan is
  fully implemented**. Treat it as a spec: if an example in it doesn't run at the end,
  the implementation is wrong (or you must flag the discrepancy — do not silently rewrite
  the README's promises). Regenerate `README.md` from it as the final step
  (`devtools::build_readme()`).

## 1. Review verdict (context for everything below)

The engine (task grid, per-task L'Ecuyer streams, checkpoint/resume, retention, model bank,
mirai transport) is solid and well-tested. The problems are one layer up:

1. **Contract inconsistencies** between the two built-in fitters and the documented interface
   (matrix orientations) that will silently corrupt any third fitter (Workstream A).
2. **The scientific vocabulary is wrong for the target audience.** "bias" and "rmse" are
   *prediction* metrics computed against the observed response (on the *training* set when no
   test set exists), while every methodologist reading this package expects the Morris,
   White & Crowther (2019) estimator-performance measures (bias of the point estimate,
   empirical SE, coverage — with MCSE). Truths are never recorded in the summary, so
   parameter-recovery analysis and `plot_recovery()` are impossible for generator-drawn
   truths (Workstream E).
3. **Interop landmines:** exported S7 generics named `fit`, `compute`, `log_lik`,
   `diagnostics` mask `generics::fit` (tidymodels), `dplyr::compute`, and
   `rstantools/brms::log_lik` for anyone who loads bayesim next to the packages this audience
   uses daily (Workstream B).
4. **Teaching material is built on a fake fitter.** Every executable example uses
   `MockFitter`, which ignores the model and fabricates draws. For a package whose stated
   purpose is education in principled Bayesian workflows, the executable path must be real
   inference (Workstream D).
5. Transport should move to purrr's mirai integration (`purrr::map()` + `in_parallel()`),
   which also forces a healthy simplification of cross-daemon error handling (Workstream C).

Work in the order A → B → C → D → E → F → G. H (cleanup) can be interleaved.

---

## Workstream A — Fix the fitter contract (CRITICAL, verified defects)

### A1. log_lik orientation is contradictory across the codebase

Verified state:
- `fit.R` generic doc (`R/fitter.R:183`): "S x N ... draws as rows ... brms/loo convention". Correct target.
- `BrmsFitter` `log_lik` returns `brms::log_lik()` = S x N. Correct.
- `MockFitter` `log_lik` (`R/fitter.R:557`) returns `matrix(nrow = n_obs, ncol = n_draws)` = **N x S. Wrong.**
- `validate_fitter()` smoke test (`R/fitter.R:965`) asserts `nrow(ll) == n` (observations as
  rows) — it **enforces the wrong orientation** and would reject a correct fitter whenever
  n_draws != n_obs.
- `build_loo_context()` doc comment (`R/worker.R:404`) says "N x S, observations in rows" while
  `relative_eff_from_chains()` (`R/worker.R:454`) correctly documents and requires S x N.
- `ElpdTestMetric` (`R/metrics-extended.R:604`) sums per-*column* log-mean-exp — correct only
  for S x N, i.e. correct for brms, silently wrong for MockFitter.

Fix: **S x N (draws x observations) everywhere** — the posterior/loo convention.
- Fix MockFitter's `log_lik`.
- Fix `validate_fitter()` to assert `ncol(ll) == n_obs` and, when the fitter reports draw
  count, `nrow(ll) == n_draws`. Make the smoke-test data deliberately non-square
  (n_obs != n_draws) so orientation bugs cannot pass.
- Fix the stray comment in `build_loo_context()`.
- Test: a deliberately transposed fitter fails `validate_fitter(smoke_test = TRUE)`;
  `elpd_test_metric` on MockFitter with n_obs != n_draws returns `n_obs` correctly.

### A2. predict_fit()$predicted_samples orientation disagrees between fitters

- Generic doc (`R/fitter.R:129`): `predicted_samples`: "N x S".
- MockFitter returns N x S; **BrmsFitter returns `posterior_predict()` output = S x N**
  (`R/brms-fitter.R:441`).
No built-in metric touches `predicted_samples` today, which is why nothing failed — it is a
landmine for user metrics. Fix: standardize on **S x N** here too (same convention as
log_lik/epred; document prominently in the Fitter contract), fix MockFitter, extend
`validate_fitter()` with a non-square orientation assertion.

### A3. BrmsFitter diagnostics ignore everything except fixed effects

`extract_brms_diagnostics()` (`R/brms-fitter.R:504`) computes rhat/ESS from `summary$fixed`
only. Hierarchical models with non-converged group-level SDs/correlations or `sigma` report
clean diagnostics. Fix: compute over **all** parameters via
`posterior::summarise_draws(as_draws_array(fit), rhat, ess_bulk, ess_tail)` (exclude `lp__`);
keep the same output field names (`rhat_max`, `ess_bulk_min`, `ess_tail_min`). Also replace the
brittle `get("control_args", asNamespace("brms"))` reach-in for max_treedepth with a documented
fallback (read `fit$fit@stan_args[[1]]$control$max_treedepth` via tryCatch, default 10L).
Test: a varying-intercept fit with deliberately few iterations reports rhat_max > threshold
even when the fixed effects look fine (skip_on_cran, cmdstan-gated).

### A4. loo_fit(BrmsFitter) drops r_eff

`R/brms-fitter.R:479` calls `loo::loo(ll)` without `r_eff`, inconsistent with the careful
chain-aware PSIS in `build_loo_context()`, giving slightly wrong p_loo/SEs vs `brms::loo()`.
Reuse `relative_eff_from_chains()` here. Test: parity of `elpd_loo` fields with
`brms::loo(fit)` on the fixture fit (tolerance 1e-6), extending
`tests/testthat/test-loo-metrics-parity.R`.

## Workstream B — API surface and naming (breaking, do once, do now)

### B1. Rename colliding generics

Renames (S7 generic + all methods + docs + vignettes + tests):

| current | new | why |
|---|---|---|
| `fit()` | `fit_model()` | masks `generics::fit` → breaks tidymodels users |
| `compute()` | `compute_metric()` | masks `dplyr::compute` |
| `log_lik()` | `log_lik_matrix()` | masks `rstantools`/`brms::log_lik` — calling `log_lik(brmsfit)` errors after `library(bayesim)` |
| `diagnostics()` | `fit_diagnostics()` | generic noun, collision-prone |
| `predict_fit()`, `predict_epred()`, `extract_draws()`, `loo_fit()` | keep | no realistic collisions |

Do not keep aliases. Test: `library(brms); library(bayesim); log_lik(some_brmsfit)` works in a
test (i.e. bayesim no longer exports `log_lik` at all).

### B2. Contract the error API

Eight error constructors + three predicates are exported. Extension authors need to *raise*
data/fit/metric errors and *classify* caught ones. Keep exported:
`bayesim_data_error`, `bayesim_fit_error`, `bayesim_metric_error`, `bayesim_config_error`,
`is_bayesim_error`, `is_fatal_error`, `is_recoverable_error`. Make internal:
`bayesim_error` (base constructor), `bayesim_contract_error`, `bayesim_checkpoint_error`,
`bayesim_internal_error`. Update `_pkgdown.yml` and the extension vignette.

### B3. Consolidate duplicate fitter validation

`validate_fitter()` (`R/fitter.R:662`, internal, has the smoke test) and
`validate_fitter_interface()` (`R/contracts.R:467`, exported, = `check_fitter_class`) overlap.
Merge into **one exported** `validate_fitter(fitter, smoke_test = FALSE)` (absorbing the class
checks), delete `validate_fitter_interface` and `check_fitter_class`, update
`validate_simulation_config()` callers, `_pkgdown.yml`, and the custom-fitters vignette. Same
treatment for `validate_metric_interface` → exported `validate_metric(metric)`.

### B4. Simplify the config knobs

- Merge `chunk_size` into `checkpoint_every` (they are `min()`-ed into one batch size anyway,
  `R/simulate.R:230`). One knob: `checkpoint_every`. Remove `chunk_size` and the deprecated
  `max_in_memory` entirely (pre-release; no shim).
- Change the data-generator signature from `(data_spec, seed, task_ctx)` to
  `(data_spec, task_ctx)`. The `seed` argument is documented as a footgun ("retained ... but
  NOT used to re-seed", `R/generators.R:11`) and `task_ctx$seed` already exists for backends
  that need an integer. Update: the three factories, `generate_template_data()`,
  `simulation_config()` validation, all vignettes/README examples, tests.
- Exclude `retain`, `max_errors`, and `checkpoint_format` from `as_config_spec()` /
  the config fingerprint (`R/simulation-config.R:435`): they are runtime policy, not
  simulation identity — changing retention must not invalidate resume.
  Test: two configs differing only in `retain` have equal fingerprints; resume across a
  `retain` change works.

### B5. S7 property validators vs constructor validation

`simulation_config()` re-implements every S7 property validator by hand
(`R/simulation-config.R:192-344`). Keep the friendly cli errors in the constructor but delete
the duplicated logic where the S7 validator already covers it (single source of truth per
rule). Low priority; do opportunistically while in the file for B4.

## Workstream C — purrr + mirai transport

Adopt purrr's mirai integration (purrr >= 1.1.0, `in_parallel()`, carrier-crated functions)
as the dispatch layer. mirai remains the daemon engine (`mirai::daemons()`,
`mirai::everywhere()` for the model bank and `daemon_setup` stay exactly as they are).

### C1. Make `run_task_safe()` total, then swap the transport

1. Change `run_task_safe()` so it **never throws**: catch fatal conditions too and return a
   failed task result whose `error` carries `fatal = TRUE` and the full condition class chain
   (`class(e)`). The controller (`execute_tasks()`) re-raises a reconstructed condition after
   collecting a batch when any result is fatal. This removes the entire cross-boundary
   condition-restoration machinery: `restore_mirai_condition()` and the
   miraiError/errorValue inspection loop in `run_batch()` (`R/simulate.R:377-500`) are deleted.
   Transport errors (daemon death) are whatever purrr/mirai raise — let them propagate as
   fatal; wrap with a clear message.
2. Replace `run_batch()`'s `mirai::mirai_map()` branch and the sequential `lapply()` branch
   with a **single** code path:

   ```r
   purrr::map(
     batch_tasks,
     purrr::in_parallel(
       \(task, run_task_safe, config_spec, fitter, metrics, retain) {
         run_task_safe(task, config_spec, fitter, metrics, retain)
       },
       run_task_safe = run_task_safe,
       config_spec = config_spec, fitter = fitter,
       metrics = metrics, retain = retain
     ),
     .progress = progress
   )
   ```

   Notes for correctness (verified against purrr docs 2026-07):
   - `in_parallel()` crates via `carrier::crate()`: the lambda's environment is stripped, so
     `run_task_safe` **must be passed as a crated dependency** (as above), not referenced
     bare. Its enclosing environment is the bayesim namespace, which resolves on daemons
     because bayesim is installed there (same requirement as today).
   - With no daemons set, purrr **falls back to sequential automatically** — delete the
     `daemons_set()` branch in `run_batch()` (keep the check in `execute_tasks()` for the
     model-bank/`daemon_setup` `everywhere()` ship).
   - Add `purrr (>= 1.1.0)` and `carrier` to Imports; keep `mirai` in Imports.
3. Unify progress: today there are two competing systems (a cli bar updated per batch in
   `execute_tasks()` + `m[mirai::.progress]`). Keep exactly one: purrr's `.progress` (named,
   e.g. `.progress = "bayesim tasks"`) inside batches; drop the outer cli bar or reduce it to
   per-batch checkpoint messages. Decide and implement one; do not ship both.
4. Tests: existing determinism tests (sequential == daemons(2)) must pass unchanged on the new
   transport; fatal-error propagation test (a `bayesim_config_error` raised inside a task
   under daemons stops the run with that condition class); the "one bank ship per run" F6
   test must still pass.

### C2. `workers` convenience argument

`run_simulation(config, workers = NULL)`: when non-NULL, call `mirai::daemons(workers)`
before task execution and `mirai::daemons(0)` in `on.exit()` — but **only when no daemons
were already set** (respect user-managed daemons; error if both `workers` and existing
daemons are present). This must happen *before* the model-bank `everywhere()` ship. Document
that `workers` is the simple path and `mirai::daemons()` the advanced/HPC path (remote
daemons, TLS, etc.). Test: `run_simulation(config, workers = 2)` matches the sequential
summary and leaves `mirai::daemons_set()` FALSE afterwards.

## Workstream D — Fitters: real inference for teaching, raw Stan for research

### D1. `LinearRegressionFitter()` — the executable-docs fitter (new)

An exact conjugate Normal–Inverse-Gamma Bayesian linear regression fitter:
`y ~ N(X beta, sigma^2)`, NIG prior (user-settable `prior_mean`, `prior_precision` scalar,
`a0`, `b0`; sensible weak defaults), formula taken from `fit_spec$formula` (base R formula,
default `response ~ .`), i.i.d. joint posterior draws of `(beta, sigma)` via the standard
NIG posterior (draw `sigma^2` from Inv-Gamma, then `beta | sigma^2` from the multivariate
normal; `n_draws` property, default 1000L). Milliseconds per fit, zero Stan, **real
posterior**. Implements the full contract including `log_lik_matrix` (exact normal density),
`predict_epred` (X %*% t(beta_draws), S x N), `predict_fit` (adds noise), `loo_fit`,
`fit_diagnostics` (rhat 1, ESS = n_draws — i.i.d. draws; divergent 0).
- Draw names must follow the plain-name convention (`Intercept`, `<coef>`, `sigma`) so
  `resolve_draw_columns()` and the generators' cleaned names work out of the box.
- This fitter replaces MockFitter in **every** user-facing example (README, all vignettes,
  pkgdown snippets). It is the package's teaching backbone: students see real posteriors,
  real coverage ≈ 0.95, real SBC uniformity, instantly.
- MockFitter: demote to internal (`@keywords internal`, un-export; keep for engine tests).
  Its "fit_idx > 1" warning hack disappears from user space with it.
- Tests: analytic posterior mean/cov vs closed form (tolerance); coverage_metric over 200
  replicates ≈ 0.95 ± MCSE; SBC ranks uniform (chi-square, generous alpha) with a
  fixed-truth-from-prior generator; `validate_fitter(smoke_test = TRUE)` passes.

### D2. `CmdStanFitter()` — user-supplied Stan files (new; user-requested)

For researchers who want their own Stan programs without brms:

```r
CmdStanFitter(
  stan_file  = "model.stan",          # or stan_code = "..."; also overridable per
                                      # fit_grid row via a `stan_file` column
  stan_data  = function(data_bundle, fit_spec) list(N = ..., y = ..., X = ...),
  log_lik    = "log_lik",             # name of the generated-quantities log-lik vector, or NULL
  epred      = NULL,                  # optional name of the epred/mu GQ matrix or vector
  chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L,
  adapt_delta = NULL, max_treedepth = NULL, parallel_chains = 1L,
  init = NULL, ...                    # passthrough to $sample()
)
```

Design decisions (implement as stated):
- **Compilation:** `cmdstanr::cmdstan_model()` caches compiled binaries by file hash natively
  — no model-bank integration needed. Compile lazily on first `fit_model()` call per distinct
  stan file, memoised in a session option keyed by file hash (same pattern as the model bank,
  shipped to daemons is NOT needed — each daemon compiles-or-cache-hits on first use; set
  `cmdstanr_write_stan_file_dir` to `result_path/stan_binaries` when result_path is set, as
  `build_model_bank()` already does, so daemons share the binary cache via the filesystem).
- **fit_model:** run `$sample(data = stan_data(data_bundle, fit_spec), seed = seed, refresh = 0, ...)`;
  capture cmdstanr warnings into the fit-result `warnings` (same `withCallingHandlers`
  pattern as BrmsFitter).
- **extract_draws:** `fit$draws()` → posterior draws, excluding `lp__` and the declared
  log_lik/epred GQ names from the default variable set (they are not parameters of interest).
- **log_lik_matrix:** extract the declared GQ vector as an S x N matrix; error with
  `bayesim_config_error` naming the convention if the variable is absent and a metric needs
  it. `supports_log_lik`/`supports_loo` auto-set from whether `log_lik` is NULL.
- **fit_diagnostics:** `fit$diagnostic_summary()` (divergences, treedepth, ebfmi) +
  `posterior::summarise_draws` rhat/ESS extrema, same field names as BrmsFitter.
- **newdata prediction is out of scope** (raw Stan has no newdata semantics):
  `supports_predictions = FALSE` unless `epred` is given, in which case `predict_epred`
  returns the in-sample GQ matrix and `predict_fit` is unsupported. Document this limitation
  clearly — test-set metrics require brms or a custom fitter.
- Vignette section (custom-fitters) with a complete worked example: a hand-written
  eight-schools or linear-regression `.stan` file with a `log_lik` GQ block, wired into an
  SBC study. Ship the `.stan` file in `inst/stan/`.
- Tests: cmdstan-gated; fit + draws + log_lik parity with an equivalent brms model
  (statistical, not bitwise); missing-GQ error message; works under `workers = 2`.

### D3. brms fit_grid ergonomics

Hand-building list-columns (`fit_grid$formula <- list(...)`) is the ugliest UX in the package
(see case-studies vignette). Add:

```r
brms_model(formula, family = NULL, prior = NULL, stanvars = NULL, stan_file = NULL)  # one spec
model_grid(gaussian = brms_model(y ~ x, gaussian()),
           student  = brms_model(y ~ x, student()))
# -> tibble with columns: model (names), formula/family/prior/stanvars list-columns
```

`model_grid()` also accepts `CmdStanFitter`-style specs later; for now brms only. Validate at
construction (formula coercible, family resolvable) so errors surface before compilation.
Update the case-studies and brms vignettes to use it. Tests: grid structure, bank dedup by
hash still works, name column lands in the summary as `fit_model`.

## Workstream E — Metrics and analysis: speak the language of simulation studies

This is the highest-leverage scientific change. Terminology and structure follow Morris,
White & Crowther (2019), *Using simulation studies to evaluate statistical methods* (Stat Med)
— cite it in the docs — and rsimsum for MCSE formulas.

### E1. Record the truth

Task results must carry the data-generating truth: in `run_task()`, after data generation,
add `true_params` to the task result as flattened `truth__<param>` columns (they are tiny —
always retained, like metrics). Without this, parameter-recovery analysis is impossible when
truths vary per task (prior-predictive/IFS generators), and `plot_recovery()`
(`R/analysis.R:315`) is **currently broken by design** — it looks for `truth`/
`true_params__*` columns that nothing ever writes. Test: summary of an IFS study contains
`truth__<param>` varying across tasks; `plot_recovery()` returns a ggplot without error.

### E2. Split estimation metrics from prediction metrics

- Rename the prediction-space metrics honestly: `rmse_metric` → `pred_rmse_metric`,
  `bias_metric` → `pred_bias_metric`, `mae_metric` → `pred_mae_metric`,
  `mse_metric` → `pred_mse_metric`. All four must **refuse to silently fall back to the
  training set**: default to test data, and when `data_bundle$test` is NULL return NA with a
  one-time warning naming the fix (provide a test set) — in-sample prediction error presented
  as "rmse" is a trap for exactly this package's audience.
- The parameter-recovery path is: `posterior_summary_metric()` (per-task point estimates +
  intervals) + `truth__*` columns (E1) + the new aggregation layer (E3). Per-task
  "parameter bias" does not need its own metric — bias/empSE/coverage are *across-replicate*
  performance measures, and computing them per task then averaging is exactly the mistake
  (mean of per-task RMSEs != RMSE) this design removes.

### E3. `performance_measures()` — the Morris et al. layer (new, exported)

```r
performance_measures(result, estimand = NULL, estimator = "mean", by = NULL)
```

For each parameter (default: all with `truth__*` + `posterior_summary__*` columns present)
and condition cell: **bias**, **empirical SE**, **MSE**, **coverage** (from
`coverage_metric()` or recomputed from `q_lower`/`q_upper` vs truth), **average model SE**
(mean posterior sd), and `n_sim` — each with its Morris et al. MCSE. Returns a tidy tibble
(one row per condition x parameter x measure or wide per measure; pick tidy-long with a
`measure` column and document it). This is the function a methods paper actually needs; make
it the centerpiece of the analysis docs.

### E4. Metric summary metadata replaces name heuristics

`summarize_simulation()` currently guesses aggregation by column-name regex
(`is_coverage_like()` treats *any* 0–1 column with ≤3 unique values as a proportion;
`rmse_mcse()` re-squares per-task RMSE values, producing an MCSE for a quantity
(sqrt-mean-MSE) that is not the reported `_mean` column — `R/analysis.R:127-151`). Metrics
are S7 objects; let them declare it: add a `summary_type` property to `Metric`
("mean", "proportion", "none") consumed by `summarize_simulation()`; drop both heuristics and
the special rmse branch (the mean of per-task pred-RMSEs gets a plain sd/sqrt(n) MCSE, which
is now correct for what is reported). Unknown/user columns default to "mean".

### E5. Metric library pruning and fixes

- **Delete `rstar_metric`** (`R/metrics-extended.R:390` returns NA unconditionally — a
  placebo export). Note it in NEWS; R* needs caret+ranger and per-chain draws; it goes on the
  backlog (I2) instead.
- `.default_prior_vars()` (`R/generators.R:315`) defaults `vars_of_interest` to population
  effects only. For SBC that silently excludes `sigma`/auxiliary parameters. Include
  distributional parameters (at minimum `sigma` when present in draws) and document how to
  override.
- Coverage/rank/pos_prob/posterior_summary already use `resolve_draw_columns()` — keep.
- Remove the vestigial `mean` field from `rank_metric` output (kept "for back-compat" in a
  pre-release package, `R/metrics-extended.R:317`) and the `references` bundle field
  (documented, always NA, consumed by nothing — remove from contracts, generators, docs).

### E6. `true_params` becomes optional

`validate_data_bundle()` hard-requires `true_params` + `vars_of_interest`
(`R/contracts.R:155-197`). Pure model-comparison/predictive studies on simulated-but-
truth-free or real data have no truths. Make both optional (NULL allowed, jointly);
truth-dependent metrics already degrade to NA. Keep the setequal check when present.
Test: a config with a truth-free generator + elpd_loo/pred_rmse metrics runs end-to-end.

### E7. Plots that earn their place

- `plot_coverage()` (`R/analysis.R:364`) puts a continuous `coverage_mean` on `geom_bar()` —
  effectively broken. Redesign: coverage (with MCSE error bars from E3) per condition/
  parameter, point-range plot, nominal line.
- `plot_recovery()`: fixed by E1; add posterior-interval segments by default, facet by
  condition.
- Keep `plot_rank_hist`/`plot_rank_ecdf` (the ECDF bands port is good) but let
  `plot_rank_hist` draw the expected-uniform band (n_reps/bins ± binomial band) like the SBC
  package.
- All plots: shared minimal theme, colorblind-safe palette, faceting by condition columns via
  a common `by` argument. Add vdiffr snapshot tests (Suggests) for all five plot functions.

## Workstream F — Runtime UX

- **F1. Preflight**: `preflight(config)` (exported) — prints/returns: task count, grid
  dimensions, metrics and their `needs` vs fitter capabilities (surfacing mismatch warnings
  *before* a 10-hour run), whether daemons are set, model-bank compilation count expected,
  estimated total time from timing one pilot task (optional `pilot = TRUE`), estimated
  checkpoint disk footprint. Also run it automatically at the top of `run_simulation()` in
  condensed one-line form ("240 tasks = 3 data x 2 fit x 40 reps; 2 models to compile;
  4 workers").
- **F2. Better failure surfacing**: at run end, when failures exist, print a compact
  `cli` summary (n failed by error_class, first message of each class) instead of leaving
  users to discover `result$errors`. Add `failed_tasks(result)` accessor returning the errors
  tibble joined with grid columns.
- **F3. `summary()`/`print()` polish**: `print.bayesim_simulation_result` should show the
  condition grid shape and a metrics preview, not just counts; add
  `as_tibble.bayesim_simulation_result` (returns `$summary`) so tidyverse users can pipe
  results directly.
- **F4. Seed ergonomics**: `simulation_config()` requires `seed` with no default — keep
  (explicitness is a feature for reproducibility teaching) but improve the error message to
  say *why* ("bayesim derives every task's RNG stream from this one seed").

## Workstream G — Documentation, demos, teaching (the emphasis)

Rules for all vignettes: every executable chunk uses `LinearRegressionFitter` (real
inference); Stan-dependent chunks are shown with `eval = knitr::is_true(...cmdstan gate...)`
and their **precomputed results shipped** as small rds files under `inst/extdata/` (generate
via a script in `data-raw/`, checked in) so pkgdown shows real output everywhere. No chunk
that fabricates numbers.

- **G1. README**: already written (`README.Rmd`, root). Regenerate `README.md` last.
- **G2. Vignette set** (replaces the current six; reuse content where good):
  1. `getting-started.Rmd` — rewrite on LinearRegressionFitter: define generator → config →
     run (with `workers = 2`) → `summarize_simulation` → `performance_measures` → one recovery
     plot and one coverage plot. Must read as "your first simulation study in 5 minutes,
     with real posteriors".
  2. `design-of-simulation-studies.Rmd` (rewrite of simulation-study.Rmd) — the teaching
     core: aims/estimands/methods/performance-measures (Morris et al. framing), choosing
     `n_replicates` from target MCSE (show the formula and a worked calculation),
     factorial vs one-at-a-time grids, reporting standards. Executable throughout.
  3. `sbc-and-calibration.Rmd` (**new flagship**) — what SBC is (Talts et al. 2018), rank
     statistics, why thinning (the `thin = "auto"` rationale), ECDF bands (Säilynoja et al.
     2022); executable LinearRegressionFitter SBC showing a calibrated model AND a
     deliberately miscalibrated one (wrong prior) so students see both outcomes;
     brms prior-predictive + IFS section with precomputed results.
  4. `brms-studies.Rmd` — model bank mechanics, `model_grid()`/`brms_model()`, family
     comparison case study (from case-studies.Rmd), warning-conditional retention, stan_args.
     Precomputed.
  5. `custom-fitters.Rmd` — the Fitter contract table (updated names/orientations, with an
     explicit "all matrices are draws x observations" callout), the LinearFitter walkthrough
     (kept), then **CmdStanFitter** with the shipped `.stan` example, then
     `validate_fitter()`.
  6. `custom-metrics.Rmd` (split out) — Metric contract, output schema, `summary_type`,
     `custom_metric()` if implemented (see G4), externalization of large outputs.
  7. `parallel-and-hpc.Rmd` (merge of memory-management + parallel content) — purrr/mirai
     model, `workers` vs `mirai::daemons()` (including remote daemons/TLS sketch for HPC),
     `daemon_setup` for cmdstan paths, retention profiles, checkpoint/resume, disk layout of
     `result_path`.
  8. `reproducibility.Rmd` — keep (it is good); update signatures and add a short "what the
     fingerprint covers" table reflecting B4.
- **G3. pkgdown**: reorder navbar (Get started = vignette 1; Articles grouped: Learn
  [2,3], Backends [4,5], Extend [5,6], Operations [7,8]); update reference index for every
  rename (B1–B3, D, E); fix the stale `news` release link ("Version 1.0.1 on CRAN" — the
  package is not on CRAN; remove the releases block); add a hex-logo placeholder slot
  (`man/figures/logo.png` — flag to the user that a logo is wanted, do not generate one).
- **G4. Lower the extension floor (evaluate, then implement if clean)**: a function-based
  metric helper so casual users never write S7:

  ```r
  custom_metric("my_rmse", needs = "predictions", summary_type = "mean",
                fn = function(fit_result, data_bundle, context, task_ctx) list(value = ...))
  ```

  Implement as an S7 `FunctionMetric` class wrapping the closure (closure must be
  crate-safe: document that it is shipped to daemons and must be self-contained). If the
  crating constraint makes this fragile under C1, document the S7 route better instead and
  record the decision in NEWS — do not ship a helper that breaks under `workers = 4`.
- **G5. Roxygen hygiene**: metric docs get the defining formula and the reference
  (Morris/rsimsum for performance measures; Vehtari et al. for LOO; Talts/Säilynoja for SBC);
  the Fitter/Metric contract pages become the single canonical statement of the interfaces
  (vignettes link, don't restate).

## Workstream H — Cleanup batch (single pass)

- `enrich_summary_with_grid_columns()` coerces task_grid-spec values to character
  (`R/simulate.R:705`, `as.character(value)`) — preserve atomic types (numeric stays
  numeric); NA only for non-scalar.
- Progress-bar variable `pb` is used outside its `if (progress && n_pending > 0)` guard at
  `R/simulate.R:367` (`cli_progress_done(id = pb)` errors if n_pending == 0 and progress
  TRUE... verify and fix; likely moot after C1's progress unification).
- Drop `dplyr` from Imports if the F7 replacement landed (verify; `grep -rn "dplyr::" R/`),
  drop `MASS` from Suggests if only the vignette LinearFitter uses it (replace with
  `stats::rnorm`-based multivariate draw or keep MASS and note why).
- `bayesim_example_data_generator` (`R/example-data-generator.R`) — unused in docs; either
  use it in getting-started or delete (flag to user first).
- `.Rbuildignore`: ensure `^PLAN.*\.md$` covers the new plan files (already updated by the
  reviewer — verify).
- NEWS.md: rewrite the 2.0.0 entry to reflect the post-plan reality (renames, new fitters,
  performance_measures, purrr transport). One coherent entry, not a changelog of the review.

## Workstream I — Future work (backlog; do NOT implement, record in a `docs/ROADMAP.md`)

- I1. Backend-agnostic prior-predictive/IFS generators (currently brmsfit-only).
- I2. R* diagnostic metric via `posterior::rstar` (needs per-chain draws — unlocked once
  extract_draws returns a posterior draws object; consider migrating the draws contract to
  `posterior::draws_df` with `.chain` in 2.2, which also removes
  `relative_eff_from_chains()`'s fit-object reach-in).
- I3. Adaptive/sequential designs (stop conditions on MCSE targets).
- I4. Quarto report template: `report(result)` renders a standard study report
  (design table, performance measures, SBC panels).
- I5. Study-level power/precision helper (invert the MCSE formula: replicates needed for a
  target coverage MCSE).
- I6. HPC recipes: SLURM + `mirai::daemons(url = ...)` worked example.
- I7. `targets` integration vignette (bayesim run as a targets pipeline step).
- I8. Optional parquet (nanoparquet) checkpoint/summary format for very large studies.
- I9. bayesfam integration vignette once `PLAN-bayesfam.md` lands (SBC as the acceptance
  test for every custom family).

## Definition of done

- [ ] A1–A4 orientation/diagnostics fixes with non-square smoke tests; transposed fitter is
      rejected by `validate_fitter()`.
- [ ] B renames complete; `log_lik` no longer exported; tidymodels/dplyr masking test green.
- [ ] C: single purrr `in_parallel()` transport; determinism (seq == workers=2) and fatal-
      propagation tests green; exactly one progress system.
- [ ] D: `LinearRegressionFitter` (analytic tests + SBC uniformity) and `CmdStanFitter`
      (cmdstan-gated tests incl. brms parity) exported and documented; MockFitter internal;
      `model_grid()`/`brms_model()` shipped.
- [ ] E: `truth__*` columns recorded; `performance_measures()` with MCSE; pred_* renames;
      `summary_type` metadata; rstar deleted; true_params optional.
- [ ] F: `preflight()`, failure summary, print/as_tibble polish.
- [ ] G: eight vignettes as specified, precomputed Stan results in `inst/extdata/`, pkgdown
      restructured, `README.md` regenerated from `README.Rmd` and its examples verified to run.
- [ ] `devtools::check()`: 0 errors / 0 warnings / 0 notes. SBC acceptance test green under
      `workers = 2`.
- [ ] NEWS.md 2.0.0 entry rewritten; ROADMAP.md created; deletion candidates flagged to the
      user, nothing deleted without confirmation.
