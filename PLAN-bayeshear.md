# bayeshear — retirement / upstreaming plan

**Audience:** an implementing agent working in the `sims1253/bayeshear` repository.
**Context:** bayeshear was a metrics companion package for bayesim 0.x/1.x. bayesim 2.0 is a
ground-up rewrite that no longer depends on bayeshear, and the metrics bayeshear provided are
now available upstream (loo, posterior, brms) or reimplemented correctly inside bayesim.
The working hypothesis (from the package owner) is that **bayeshear is no longer useful and
should be retired**. This plan is an audit-then-retire brief: verify the hypothesis function
by function, upstream the few pieces worth keeping, then archive.

**Hard constraint carried over from bayesim's review:** bayeshear's LOO-RMSE / LOO-R²
implementations were found to be **statistically wrong** (invented elpd-to-RMSE identities
that hold only for fixed-variance Gaussians). They must never be ported anywhere. The
correct constructions live in bayesim 2.0 (`loo::E_loo()`-based, matching
`brms::loo_R2()` / `loo_predict()`) with parity tests against brms.

## Process requirements

- Never run state-modifying git commands; the user handles git. Read-only git is fine.
- Never delete files without explicit user confirmation — this plan *ends* in archival, but
  the archival step (repo settings, final commit) is the user's.
- Produce the audit table (step 1) as `AUDIT.md` in the repo root before touching code.

## Step 1 — Function-by-function audit

For every exported function in bayeshear, record in `AUDIT.md`: name, what it computes, and
its disposition with evidence. Expected dispositions:

| disposition | criterion | expected examples |
|---|---|---|
| **superseded-upstream** | an actively maintained equivalent exists in loo / posterior / brms | LOO predictive metrics → `loo::loo_predictive_metric()` (mae/rmse/sq. error via E_loo); R² → `brms::loo_R2()` / `brms::bayes_R2()`; R* → `posterior::rstar()`; rank-normalization / u-scale → `posterior::u_scale()` |
| **superseded-bayesim** | reimplemented (correctly) in bayesim 2.0 | SBC rank machinery, ECDF bands (`adjust_gamma` port), coverage/bias summaries |
| **known-wrong** | the incorrect LOO-RMSE/R² family | flag prominently; never port |
| **unique-keep** | no upstream equivalent AND still scientifically useful | candidate: the multi-chain (L > 1) `adjust_gamma_simulate()` simultaneous-band calibration that bayesim's port stubbed out (bayesim falls back to the L = 1 band, which is conservative). Verify whether `posterior`/SBC package now covers this; if not, this is the one concrete upstreaming candidate into bayesim (`R/sbc-bands.R`) |

Verify "superseded-upstream" claims against the **current installed versions** of loo,
posterior, and brms (read their docs/NEWS — do not trust memory; loo gained
`loo_predictive_metric()` and `E_loo()` improvements across 2.x).

## Step 2 — Upstream the survivors

For each `unique-keep` item, write a self-contained patch **as a diff/branch against
sims1253/bayesim** (do not modify bayesim from this repo without coordination — attach the
patch to the audit): most likely a faithful `adjust_gamma_simulate()` (Monte-Carlo gamma
calibration for L > 1 chains) plus `u_scale()` via `posterior::u_scale`, wired into
bayesim's `adjust_gamma()` L > 1 branch, with a test comparing L = 1 exact vs L > 1
simulated bands on synthetic uniform ranks.

## Step 3 — Retirement

1. Add a prominent README banner: package is superseded; point to bayesim 2.x for simulation
   metrics and to loo/posterior/brms for the individual quantities, with a short
   old-name → new-home mapping table (generated from AUDIT.md).
2. Add the same mapping to a final NEWS entry.
3. Do NOT half-maintain: no deprecation-warning release cycle is needed for a GitHub-only
   research package — one final informative state is enough.
4. Hand back to the user with: AUDIT.md, the bayesim patch (if any), and the recommendation
   to archive the GitHub repository (user action).

## Definition of done

- [ ] AUDIT.md covers 100% of exports with evidence links (upstream function + version).
- [ ] The known-wrong LOO metrics are explicitly flagged as never-port.
- [ ] Any unique functionality exists as a tested patch against bayesim, not as a reason to
      keep bayeshear alive.
- [ ] README/NEWS point every old entry point to its new home.
- [ ] Repo is ready for the user to archive.
