# Reproducibility

## Reproducibility guarantees

With the same configuration, package and backend versions, and platform,
bayesim reproduces scientific outputs across sequential, parallel, and
resumed runs. Wall-clock timings differ. Custom generators, fitters, and
metrics must use the supplied RNG state and avoid mutable external
state.

Changing platforms or backend versions can change numerical results.
Check the results again after an upgrade; the seed alone cannot
guarantee that different implementations produce equivalent inference.

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

Because the RNG is per-task and stream-derived, sequential runs and
mirai daemon runs produce the same scientific outputs and canonical task
outcomes. Volatile timing fields are the only difference, so the example
drops them before comparing:

``` r

# Sequential
seq_result <- run_simulation(
  config, resume = "never", progress = FALSE, verbose = FALSE
)

# Parallel via mirai
mirai::daemons(2)
par_result <- run_simulation(
  config, resume = "never", progress = FALSE, verbose = FALSE
)
mirai::daemons(0)

# Timing columns record wall-clock time, so drop them before comparing.
drop_timing <- function(df) df[!grepl("^timing_", names(df))]
identical(
  drop_timing(seq_result$summary),
  drop_timing(par_result$summary)
) # TRUE
#> [1] TRUE
```

#### Determinism under resume

When a run is interrupted and resumed, the task grid (including each
task’s RNG stream) is recomputed deterministically from the seed, and
only terminal task statuses are copied from the checkpoint. Resuming
therefore completes the remaining tasks with the exact RNG streams they
would have had in a full run. The resumed study then matches a single
uninterrupted run once volatile timing fields are excluded.

### The config fingerprint

`config_fingerprint(config)` returns the SHA256 study identifier stored
in the checkpoint manifest. Resume requires a matching fingerprint.

The fingerprint includes the data, fit, and explicit task grids;
replicate count; seed; generator signature; and fitter and metric
classes, properties, and available package versions. A generator
signature includes argument names and a body hash. For a package
function it also records the package, function name, and version when an
unambiguous reference is available. For a closure it hashes referenced
values bound directly in its local environment.

Editing a generator body or a captured local value therefore changes the
fingerprint. Global variables, inherited environment bindings, and
external file contents are not tracked by that signature. Keep study
inputs in the grids or explicit local bindings, and preserve external
inputs separately. The fingerprint is a compatibility check, not a
complete dependency record.

Runtime settings are excluded: output path, retention, checkpoint
settings, error limits, adaptive stopping (`stop_on`), and daemon setup.
You can adjust these without changing the study identifier. Resume still
checks retention compatibility; it cannot recover artifacts that an
earlier run discarded.

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
- Upgrading Stan/brms/cmdstanr, which can change numerical behavior.
