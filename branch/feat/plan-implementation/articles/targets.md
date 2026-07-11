# Running bayesim with targets

``` r

library(bayesim)
```

## Introduction

[targets](https://docs.ropensci.org/targets/) is a Make-like pipeline
engine for R: it caches the return value of each step and only re-runs a
step when its inputs change. That maps onto bayesim naturally because
every simulation has a stable **config fingerprint** — a SHA256 hash of
the data/fit grids, the generator/fitter/metrics specs, `n_replicates`,
and the seed (see
[`vignette("reproducibility")`](https://sims1253.github.io/bayesim/articles/reproducibility.md)).
The exported
[`config_fingerprint()`](https://sims1253.github.io/bayesim/reference/config_fingerprint.md)
accessor is useful as an audit ID and as a compact dependency for
downstream targets.

`targets` is not a hard dependency of bayesim. The pipeline code below
is `eval = FALSE`; install `targets` separately if you want to run it.

## A minimal targets pipeline

A typical workflow defines one dynamically branched target per
simulation condition (so conditions rebuild independently) or one target
for the whole config. `targets` already hashes commands and dependencies
under its default `tar_cue(mode = "thorough")`; there is no
`fingerprint` argument to `tar_cue()`. The pipeline below relies on that
native invalidation and records bayesim’s fingerprint as a separate
target.

``` r

# _targets.R
library(targets)
library(bayesim)

# Build the data and fit grids once, up front.
data_grid <- data.frame(
  n = c(100, 500, 1000),
  effect = c(0.2, 0.5, 1.0)
)
fit_grid <- data.frame(model = c("baseline", "full"))

# One row per condition; we will map over it to define one target per condition.
conditions <- split(expand.grid(
  data_idx = seq_len(nrow(data_grid)),
  fit_idx = seq_len(nrow(fit_grid)),
  stringsAsFactors = FALSE
), seq_len(nrow(data_grid) * nrow(fit_grid)))

list(
  tar_target(condition, conditions, iteration = "list"),

  # One dynamic branch per (data_idx, fit_idx) condition.
  tar_target(
    sim_condition,
    {
      config <- simulation_config(
        data_grid = data_grid[condition$data_idx, , drop = FALSE],
        fit_grid = fit_grid[condition$fit_idx, , drop = FALSE],
        data_generator = my_data_generator,
        fitter = LinearRegressionFitter(n_draws = 500L),
        metrics = list(posterior_summary_metric(), coverage_metric()),
        n_replicates = 200L,
        seed = 42L
      )
      run_simulation(config, progress = FALSE)
    },
    pattern = map(condition)
  ),

  # A single-target alternative for the whole study.
  tar_target(
    config_all,
    simulation_config(
      data_grid = data_grid,
      fit_grid = fit_grid,
      data_generator = my_data_generator,
      fitter = LinearRegressionFitter(n_draws = 500L),
      metrics = list(posterior_summary_metric(), coverage_metric()),
      n_replicates = 200L,
      seed = 42L
    )
  ),
  tar_target(config_id, config_fingerprint(config_all)),
  tar_target(
    sim_all,
    {
      config_id # explicit audit dependency, also visible in the target graph
      run_simulation(config_all, progress = FALSE)
    }
  )
)
```

Notes on invalidation:

- The default thorough cue rebuilds a target when its command or
  upstream dependencies change. No custom cue is needed.
- `config_id` excludes runtime policy (`result_path`, retention,
  checkpoint settings, `max_errors`, and `stop_on`) by design, making it
  a stable study ID. `sim_all` still depends on the full `config_all`,
  so changing runtime policy can correctly affect execution and
  checkpoint behavior.

## Combining results

Once the per-condition targets have built, assemble them into a single
performance table with `tar_read()` / `tar_load()`. bayesim’s
[`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md)
expects a summary data frame, so we row-bind the per-condition summaries
first:

``` r

# In an interactive session, after tar_make():
library(targets)
# Load all per-condition results.
conditions_results <- tar_read(sim_condition)          # list, one element per branch
combined <- do.call(rbind, lapply(conditions_results, \(r) r$summary))

pm <- performance_measures(combined, estimand = "effect")
pm
```

For the whole-config target, `tar_read(sim_all)$summary` gives the full
per-task summary directly, and `summarize_simulation(tar_read(sim_all))`
produces the aggregated condition-level table with MCSEs.

## Next steps

- [`vignette("getting-started")`](https://sims1253.github.io/bayesim/articles/getting-started.md)
  for
  [`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
  and
  [`performance_measures()`](https://sims1253.github.io/bayesim/reference/performance_measures.md).
- [`vignette("reproducibility")`](https://sims1253.github.io/bayesim/articles/reproducibility.md)
  for what exactly feeds the config fingerprint (and what does not).
- [`vignette("parallel-and-hpc")`](https://sims1253.github.io/bayesim/articles/parallel-and-hpc.md)
  for the remote-daemon / SLURM side of large runs, which composes
  naturally with a `targets` driver.
