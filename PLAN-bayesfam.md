# bayesfam — hardening plan as the family library for bayesim studies

**Audience:** an implementing agent working in the `sims1253/bayesfam` repository.
**Context:** bayesfam provides custom brms response families (exotic likelihoods beyond
brms's built-ins). bayesim 2.0 (`sims1253/bayesim`) is a simulation framework whose
generators (`prior_predictive_generator()`, `ifs_generator()`) forward-sample data through
`brms::posterior_predict()`, and whose metrics consume `brms::log_lik()` and
`brms::posterior_epred()`. That makes bayesfam's role sharp: **every family bayesfam ships
must round-trip through posterior_predict / log_lik / posterior_epred correctly**, because
bayesim will use those three functions to generate data from the family, refit, and check
calibration. The synergy also gives bayesfam its acceptance test for free: a family whose
SBC ranks are uniform under bayesim has a self-consistent (rng, lpdf, mean) triple.

**Do not trust memorized brms internals** — custom-family hooks (`posterior_predict_<fam>`,
`log_lik_<fam>`, `posterior_epred_<fam>`, `custom_family()` signatures) have changed across
brms 2.x; verify against the installed brms version's custom-families vignette first.

## Process requirements

- Never run state-modifying git commands; the user handles git.
- Never delete files without explicit user confirmation.
- bayesim is installed and cmdstan is available; long tests gate on cmdstan and
  `skip_on_cran()`.

## Step 1 — Family capability audit

For every exported family, record in `AUDIT.md`: family name, link(s), and whether each hook
exists and is *tested*: `log_lik_<fam>` (pointwise, draws x obs), `posterior_predict_<fam>`
(rng), `posterior_epred_<fam>` (analytic mean, no noise). Missing epred or rng hooks make the
family unusable for bayesim's IFS/prior-predictive generators and LOO metrics — that is the
gap list this plan closes.

## Step 2 — Per-family correctness tests (the core work)

For each family, add:

1. **rng ↔ lpdf self-consistency:** large-sample draws from the rng match the density
   (Kolmogorov–Smirnov or moment checks against numeric integration of the lpdf over a
   parameter grid including boundary-adjacent values).
2. **epred correctness:** `posterior_epred` equals the analytic mean; verify numerically
   (quadrature over the lpdf) where a closed form is fragile.
3. **SBC acceptance via bayesim (the flagship test):** one `tests/testthat/test-sbc-<fam>.R`
   per family — `bayesim::prior_predictive_generator()` on a prior-only intercept model of
   the family → `BrmsFitter` (small iter, model bank on) → `rank_metric()` +
   `coverage_metric()`, ~50–100 replicates, `mirai` daemons optional. Assert rank uniformity
   (chi-square, generous alpha, fixed seed) and coverage in a wide band. This single test
   catches sign errors, wrong Jacobians, mismatched parameterizations between the Stan lpdf
   and the R-side hooks — the classic custom-family failure modes.
   Check bayesim's current API before writing these (the 2.1 plan renames `fit` →
   `fit_model` etc.; the generator/metric constructors above are stable names).
4. **Parameterization documentation:** every family's doc page states the density formula,
   the parameter meanings, link defaults, and the mean function — these are teaching
   packages now; a student must be able to check the math.

## Step 3 — CI and compatibility

- CI matrix against brms release and brms dev (custom-family breakage historically comes
  from brms dev changes); cmdstan cached.
- A single "smoke everything" job: compile one model per family (`chains = 0`) so Stan-code
  regressions surface fast without sampling.

## Step 4 — bayesim-facing deliverable

- A `families()` catalog function returning a tibble (name, support, links, has_epred,
  has_rng, sbc_status) that bayesim's future integration vignette (bayesim ROADMAP item I9)
  can render.
- One worked example in the README: a bayesim family-recovery study (does fitting family X
  to data generated from family Y bias the coefficients?) — this is the package's selling
  point for simulation studies.

## Definition of done

- [ ] AUDIT.md: full hook coverage table.
- [ ] Every family has rng/lpdf/epred tests; gaps either fixed or the family is flagged to
      the user as a removal candidate (do not delete unilaterally).
- [ ] SBC-via-bayesim acceptance tests pass for all retained families.
- [ ] CI matrix green on brms release + dev.
- [ ] `families()` catalog + README study example.
