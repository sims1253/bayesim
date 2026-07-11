# Parallel execution, HPC, and memory

``` r

library(bayesim)
```

This vignette covers running studies beyond a single core and a single
sitting: the purrr/mirai execution model, `workers` vs user-managed
[`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
(including remote/TLS daemons on a cluster), the `daemon_setup` hook,
retention profiles and batching for flat memory, and the on-disk layout
of `result_path` for checkpoint/resume.

## The execution model

bayesim dispatches tasks with
[`purrr::map()`](https://purrr.tidyverse.org/reference/map.html) +
[`purrr::in_parallel()`](https://purrr.tidyverse.org/reference/in_parallel.html);
[mirai](https://mirai.r-lib.org/) is the daemon engine underneath. There
is a single code path: with no daemons set, purrr runs tasks
sequentially in the current process; with daemons set, the same tasks
are crated and shipped to the workers. Determinism is unaffected — every
task carries its own precomputed L’Ecuyer-CMRG RNG stream, so sequential
and parallel runs produce identical results (see
[`vignette("reproducibility")`](https://sims1253.github.io/bayesim/articles/reproducibility.md)).

Task lambdas are shipped to daemons via
[`carrier::crate()`](https://rdrr.io/pkg/carrier/man/crate.html), which
strips their environment. Your data generator and any custom
fitters/metrics are crated into the transport, so they must be
self-contained: reference package functions by namespace and pass data
through `data_spec`, not through free variables captured from your
session.

## Local parallelism: `workers`

For local multi-core work, pass `workers` to
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md):

``` r

result <- run_simulation(config, workers = 4)
# mirai::daemons(4) is called for you and torn down on exit.
```

This is the recommended path for development and moderate task counts.
`workers` errors if daemons are already set (it will not silently stack
a second daemon pool on top of yours); everything below is for daemons
the `workers` argument cannot manage itself — remote nodes, TLS,
per-daemon setup.

## Remote daemons (mirai URL scheme)

mirai’s URL scheme splits the work into two roles:

- A **dispatcher** (the controller process) that calls
  `mirai::daemons(url, n)` to advertise `n` worker slots at a URL.
- **Daemons** (the workers), launched separately, that each call
  `mirai::daemon(url, ...)` to connect back and pull tasks.

This split maps cleanly onto a batch scheduler: the dispatcher runs on
the head node, the daemons on the compute nodes.

``` r

library(bayesim)
library(mirai)

# On the dispatcher: advertise 16 worker slots at a TLS-secured URL.
url <- "tls://10.0.0.1:5555"
daemons(url = url, n = 16)

# Run bayesim with workers = NULL so it uses the daemons you set up.
result <- run_simulation(config, workers = NULL)

daemons(0)
```

### A SLURM submission script

A typical pattern is one SLURM job that launches `N` daemon workers
across nodes, then runs the dispatcher R script on the head node:

``` bash
#!/bin/bash
#SBATCH --job-name=bayesim-sim
#SBATCH --nodes=4
#SBATCH --tasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=sim-%j.out

module load R
module load gcc

export BAYESIM_TLS=/scratch/$USER/bayesim-tls

DISPATCHER_URL="tls://$(hostname -I | awk '{print $1}'):5555"
N_DAEMONS=$((SLURM_NTASKS - 1))   # reserve one task for the dispatcher

srun --ntasks=$N_DAEMONS --overlap Rscript -e "
  mirai::daemon('$DISPATCHER_URL',
                tls = mirai::tls_config(
                  ca = '$BAYESIM_TLS/ca.pem',
                  key = '$BAYESIM_TLS/server-key.pem',
                  cert = '$BAYESIM_TLS/server.pem'))
" &

sleep 10  # give the daemons a moment to connect

Rscript dispatcher.R "$DISPATCHER_URL" "$N_DAEMONS"

wait
```

And the dispatcher R script:

``` r

# dispatcher.R
args <- commandArgs(trailingOnly = TRUE)
url <- args[1]
n <- as.integer(args[2])

library(bayesim)
mirai::daemons(url = url, n = n)

config <- simulation_config(
  data_grid = my_data_grid,
  fit_grid = my_fit_grid,
  data_generator = my_data_generator,
  fitter = BrmsFitter(),
  metrics = list(posterior_summary_metric(), coverage_metric()),
  n_replicates = 500L,
  seed = 42L,
  result_path = "sim-checkpoints"
)

result <- run_simulation(config, workers = NULL)

mirai::daemons(0)
```

The exact `srun`/socket plumbing depends on your cluster’s network
topology; the load-bearing idea is the dispatcher/daemon split.

### TLS configuration

Remote daemons should connect over TLS. mirai uses standard PEM
material; see the [mirai security
documentation](https://mirai.r-lib.org/articles/tls.html) for the full
reference. A self-signed CA plus server certs is enough to get started:

``` bash
# Run once, on a trusted host. Keep the private keys private.
openssl req -x509 -newkey rsa:2048 -nodes -days 365 \
  -keyout ca-key.pem -out ca.pem \
  -subj "/CN=bayesim-CA"

openssl req -newkey rsa:2048 -nodes \
  -keyout server-key.pem -out server-csr.pem \
  -subj "/CN=$(hostname)"

openssl x509 -req -in server-csr.pem -CA ca.pem -CAkey ca-key.pem \
  -CAcreateserial -out server.pem -days 365

# Distribute ca.pem, server.pem, and server-key.pem to every node
# that runs a daemon or the dispatcher.
```

Then pass the material into both sides via `mirai::tls_config()` as in
the SLURM script above.

## The daemon_setup hook

Daemons start from a fresh R process. Anything the workers need at run
time — the cmdstan install path, environment options, a warm cache —
must be set up on each daemon before tasks start. The `daemon_setup`
argument of
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
runs once per daemon via
[`mirai::everywhere()`](https://mirai.r-lib.org/reference/everywhere.html)
before the task grid ships:

``` r

config <- simulation_config(
  ...,
  daemon_setup = function() {
    # Point cmdstanr at the cluster's cmdstan install.
    cmdstanr::set_cmdstan_path("/opt/cmdstan-2.35.0")
  }
)
```

`daemon_setup` is ignored when no daemons are set, so the same config
works for local and HPC runs. For `BrmsFitter`, the model bank (compiled
binaries) is also shipped to each daemon once per run — daemons must be
set *before*
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
is called.

## Keeping memory flat

Large studies (thousands of tasks, each with a full posterior) can
exhaust memory. Two knobs keep it bounded:

### Retention profiles

The `retain` argument of
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)
selects what each task result keeps:

| value | keeps | drops |
|----|----|----|
| `"minimal"` | computed metric values | diagnostics and heavy artifacts |
| `c("metrics", "diagnostics")` (default) | metrics and diagnostics | warnings and heavy artifacts |
| `"standard"` | metrics, diagnostics, warnings | fit object, draws, predictions, data |
| `"debug"` | every retention option | nothing |
| an explicit vector | exactly the requested fields (plus metrics) | every unrequested field |

For a metrics-only sweep over thousands of conditions,
`retain = "metrics"` keeps memory flat regardless of posterior size. Use
`"debug"` only when you need all per-task artifacts, and consider the
conditional form (`retain = list(success = ..., warning = ...)`) to keep
heavy artifacts only for tasks that emitted warnings — see
[`vignette("brms-studies")`](https://sims1253.github.io/bayesim/articles/brms-studies.md).

### Batching: `checkpoint_every`

`checkpoint_every` is the single batching knob: the engine processes
tasks in batches of this size, writes a checkpoint after each batch
(when `result_path` is set), and lightens completed results per the
retention profile before starting the next batch.

``` r

config <- simulation_config(
  ...,
  result_path = "results/big_study",
  checkpoint_every = 50L,
  keep_checkpoints = 2L,
  retain = "metrics"
)
```

When using daemons, each daemon holds its own copies of in-flight task
artifacts; `checkpoint_every` bounds the in-flight work and the
retention profile bounds per-task size.

Each checkpoint is a complete cumulative snapshot so it can be validated
and resumed independently. `keep_checkpoints = 2L` (the default) prunes
older snapshots after a successful write, bounding checkpoint disk usage
while preserving one fallback if the newest snapshot is corrupted.
Increase it only when you need a longer audit history.

## Checkpoint/resume and the result_path layout

Point `result_path` at a directory (on shared scratch, for clusters) and
pass `resume = "auto"` so an interrupted run picks up where it left off:

``` r

result <- run_simulation(config, resume = "auto", workers = NULL)

# Or force a resume from a specific directory:
result <- resume_simulation("/scratch/user/bayesim-sim")
```

The directory contains everything needed to resume and audit a run:

    results/my-study/
    ├── run_manifest.json      # config fingerprint, schema version, spec summary
    ├── latest.json            # pointer to the newest complete checkpoint
    ├── checkpoints/           # atomic rds checkpoints (task grid + results)
    ├── artifacts/             # externalized large metric outputs
    │   └── metrics/
    ├── stan_binaries/         # shared compiled-model cache (Stan backends)
    └── summary.parquet        # optional sidecar (summary_format = "parquet")

The config fingerprint is written into the manifest, so a resumed run is
rejected if you change the study design (grids, seed, generator, fitter
or metric specs) between submissions. Runtime policy — `retain`,
`max_errors`, `checkpoint_every`, `keep_checkpoints`,
`checkpoint_format` — is deliberately *excluded* from the fingerprint:
you can change retention or error tolerance and still resume.

## Next steps

- [`vignette("reproducibility")`](https://sims1253.github.io/bayesim/articles/reproducibility.md)
  — why resuming and parallelizing produce identical results.
- [`vignette("brms-studies")`](https://sims1253.github.io/bayesim/articles/brms-studies.md)
  — the model bank and warning-conditional retention.
- [`vignette("targets")`](https://sims1253.github.io/bayesim/articles/targets.md)
  — running bayesim as a step in a `targets` pipeline.
