# bayesim roadmap (backlog)

Items recorded during the 2.0.0 redesign (see `PLAN.md`, Workstream I). Most
were implemented in the deferred-work pass; the remaining open items are marked
**OPEN**.

## Implemented (deferred-work pass)

- **I1.** Backend-agnostic `prior_draws_generator()` / `forward_sim_generator()`
  — work via the S7 `extract_draws`/`predict_fit` generics. NOTE: for non-brms
  fitters the "prior" is approximate (posterior-as-prior); brms users keep the
  brms-specific `prior_predictive_generator()` for full prior predictive
  coverage. A true generic prior would need a `prior_draws` S7 method (below).
- **I2.** `rstar_metric()` re-added via `posterior::rstar` (per-chain draws;
  NA + warning for chain-less fitters; gated on caret/randomForest for brms).
- **I3.** Adaptive stopping via `simulation_config(stop_on = ...)` — skips
  remaining tasks once an estimand/measure MCSE falls below target.
- **I4.** `report()` renders a Quarto study report from a template.
- **I5.** `n_replicates_for_target()` inverts the MCSE formula.
- **I6.** SLURM + remote mirai daemons, TLS, daemon_setup — now part of
  `vignettes/parallel-and-hpc.Rmd` (which also absorbed memory-management).
- **I7.** `vignettes/targets.Rmd` — bayesim as a `targets` pipeline step.
- **I8.** Optional parquet summary output (`summary_format = "parquet"`,
  nanoparquet; `read_summary()`).
- **I9.** `vignettes/bayesfam.Rmd` stub (SBC-as-acceptance-test; pending bayesfam release).

## OPEN

- **A true generic `prior_draws` S7 method.** `prior_draws_generator()` for
  non-brms fitters uses the pilot fit's posterior draws as an approximation.
  Adding `prior_draws(fitter, fit_spec, n)` as an S7 generic (with a default
  that returns the weak-prior posterior for `LinearRegressionFitter` and a
  `brms::prior_draws` wrapper for `BrmsFitter`) would make the prior-predictive
  path exact for all fitters.
- **R\* with `ranger`/`gbm` backends.** The current `rstar_metric` uses
  caret's default `randomForest` backend; expose `method` / `hyperparameters`
  (already properties on `RstarMetric`) in the constructor and document.
- **`bayesfam` integration** — the vignette stub points to the as-yet-unreleased
  `bayesfam` package; fill in the real workflow once it lands.
- **Full `targets` round-trip test** — the targets vignette is prose-only
  (eval=FALSE); a runnable end-to-end `targets` pipeline test would lock in the
  cache-cue contract.
- **Rank-histogram (binned) view** — `plot_rank_ecdf` has the ECDF + uniformity
  band; a complementary binned histogram with the discrete uniform envelope is
  still wanted for the SBC vignette.
