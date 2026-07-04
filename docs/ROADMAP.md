# bayesim roadmap (backlog)

Items deliberately deferred from the 2.0.0 redesign (see `PLAN.md`,
Workstream I). Not implemented; recorded here so the work is not lost.

## I1. Backend-agnostic prior-predictive / IFS generators
`prior_predictive_generator()` and `ifs_generator()` currently require a
brmsfit. Generalise to any fitter that exposes prior/posterior draws (incl.
`CmdStanFitter` and `LinearRegressionFitter`).

## I2. R\* diagnostic metric via `posterior::rstar`
Re-add `rstar_metric()` (deleted in 2.0.0 as a placebo — it returned NA
unconditionally). Needs per-chain draws (caret + ranger). Unlocked once
`extract_draws` returns a `posterior::draws_df` with `.chain` (consider
migrating the draws contract in 2.2, which also removes
`relative_eff_from_chains()`'s fit-object reach-in).

## I3. Adaptive / sequential designs
Stop conditions on MCSE targets — pause/stop a run once a coverage or bias MCSE
threshold is met.

## I4. Quarto report template
`report(result)` renders a standard study report (design table, performance
measures, SBC panels) to HTML/PDF.

## I5. Study-level power / precision helper
Invert the MCSE formula: given a target coverage MCSE, return the number of
replicates needed.

## I6. HPC recipes
SLURM + `mirai::daemons(url = ...)` worked example vignette.

## I7. `targets` integration
bayesim run as a `targets` pipeline step (cache invalidation via the config
fingerprint).

## I8. Optional parquet (nanoparquet) checkpoint / summary format
For very large studies where the rds checkpoint / wide summary becomes a
bottleneck.

## I9. bayesfam integration vignette
SBC as the acceptance test for every custom family, once
`PLAN-bayesfam.md` lands.

## Deletion candidates (flagged, not deleted — confirm with maintainer)
- `bayesim_example_data_generator` (`R/example-data-generator.R`): unused in
  docs; either wire into `getting-started` or remove.
- `references` bundle field: documented, always NA, consumed by nothing —
  candidate for removal from the data contract / generators / docs (the E5
  sweep left it optional but in place).
