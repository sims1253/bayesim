# Reproducibility

## Reproducibility guarantees

bayesim is designed so that simulation results are reproducible across
sequential, parallel, interrupted-and-resumed, and re-run executions,
given a fixed seed.

### What is guaranteed

- **Same seed + same package/backend versions + same platform** produces
  **identical** results (bit-for-bit), whether the run is sequential,
  uses mirai daemons, or is interrupted and resumed.
- **Across platforms** (or when the Stan/brms/cmdstanr version changes),
  results are **statistically equivalent** but may differ in the least
  significant bits due to floating-point ordering and backend
  differences.

This is the honest, qualified statement; bit-identical reproducibility
cannot be promised across Stan versions or operating systems.

### How determinism is achieved

#### Per-task RNG streams

At the start of a run, bayesim derives one L’Ecuyer-CMRG RNG stream per
task from the simulation seed. Each task restores its stream before its
data generator and fitter run, so the RNG state a task sees depends only
on its position in the grid – not on execution order, parallelism, or
which other tasks have completed.

Stochastic metrics receive deterministic sub-seeds derived from the task
seed and metric name. Adding or reordering metrics therefore does not
change the random draws used by another metric.
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
also restores the caller’s RNG kind and state when it returns.

``` r

library(bayesim)

config <- simulation_config(
  data_grid = data.frame(n = c(50, 100)),
  fit_grid = data.frame(model = "m"),
  data_generator = function(data_spec, task_ctx) {
    # consume the AMBIENT RNG state (restored by the worker before this call);
    # do NOT call set.seed() or withr::with_seed() here
    n <- data_spec$n
    x <- stats::rnorm(n)
    y <- x + stats::rnorm(n)
    list(
      train = data.frame(y = y, x = x), test = NULL, response = "y",
      true_params = c(slope = 1), vars_of_interest = "slope"
    )
  },
  fitter = LinearRegressionFitter(n_draws = 200L),
  metrics = list(posterior_summary_metric()),
  n_replicates = 2L,
  seed = 42L
)
```

#### Determinism across executors

Because the RNG is per-task and stream-derived, a study produces
identical summaries whether run sequentially or via mirai daemons:

``` r

# Sequential
seq_result <- run_simulation(config, resume = "never", progress = FALSE)

# Parallel via mirai
mirai::daemons(2)
par_result <- run_simulation(config, resume = "never", progress = FALSE)
mirai::daemons(0)

identical(seq_result$summary, par_result$summary) # TRUE (ignoring timing columns)
```

#### Determinism under resume

When a run is interrupted and resumed, the task grid (including each
task’s RNG stream) is recomputed deterministically from the seed, and
only terminal task statuses are copied from the checkpoint. Resuming
therefore completes the remaining tasks with the exact RNG streams they
would have had in a full run, producing identical results.

### The config fingerprint

Each simulation configuration is hashed into a stable fingerprint that
is written to the checkpoint manifest. On resume, bayesim verifies the
fingerprint matches, so a resumed run cannot silently use a different
config than the one that wrote the checkpoint.

### Generators and determinism

The factory generators
([`fixed_truth_generator()`](https://sims1253.github.io/bayesim/reference/fixed_truth_generator.md),
[`prior_predictive_generator()`](https://sims1253.github.io/bayesim/reference/prior_predictive_generator.md),
[`ifs_generator()`](https://sims1253.github.io/bayesim/reference/ifs_generator.md))
all consume the ambient RNG state. IFS and prior-predictive generators
additionally select their truth-draw parameter vector by a
**deterministic** index derived from `task_ctx$rep_idx` – never by
random sampling of draw indices – so SBC ranks are well-defined and
resume is reproducible.

### What can break reproducibility

- Changing the seed.
- Changing the data generator’s RNG consumption (e.g. adding an extra
  `rnorm` call) – this shifts all downstream draws for that task.
- Renaming a stochastic metric, because its name is part of its
  deterministic metric-specific seed.
- Reordering the data/fit grid (tasks are identified by grid position).
- Upgrading Stan/brms/cmdstanr – results stay statistically equivalent
  but not bit-identical.
