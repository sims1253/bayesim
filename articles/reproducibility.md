# Reproducibility in bayesim

## Introduction

Reproducibility is a requirement for rigorous scientific research. In
simulation studies, reproducibility means that running the same
simulation with the same inputs should produce identical results—every
time, on any machine, regardless of interruptions or execution mode.

bayesim provides deterministic reproducibility as a core guarantee. This
vignette explains:

- **Why** reproducibility matters
- **How** bayesim ensures reproducibility through careful RNG management
- **What** is and isn’t guaranteed to be deterministic
- **Best practices** for ensuring your simulations remain reproducible

### Why Reproducibility Matters

Simulation studies are experiments, and like any experiment, they must
be reproducible to be scientifically valid:

1.  **Verification**: Others (and future you) can verify your results
2.  **Debugging**: Deterministic behavior makes it possible to reproduce
    and fix bugs
3.  **Collaboration**: Team members can share and build upon each
    other’s work
4.  **Publication**: Reviewers and readers can replicate your findings
5.  **Longevity**: Results remain valid even as hardware and software
    evolve

Without reproducibility, simulation results are unverifiable claims.

## How bayesim Ensures Reproducibility

bayesim achieves reproducibility through four mechanisms:

### 1. L’Ecuyer-CMRG RNG for Parallel-Safe Streams

bayesim uses the **L’Ecuyer-CMRG** (Combined Multiple Recursive
Generator) random number generator instead of R’s default
Mersenne-Twister. This choice is deliberate:

- **Multiple independent streams**: L’Ecuyer-CMRG can generate multiple
  independent RNG streams from a single seed
- **Reproducible parallel execution**: Each task gets its own stream,
  ensuring results do not depend on execution order
- **Long period**: The generator has a long period, preventing overlap
  between streams

``` r
library(bayesim)

# bayesim automatically sets L'Ecuyer-CMRG when you run a simulation
# You can also set it manually for custom code:
setup_global_rng(42)

# Check the RNG kind
RNGkind()
#> [1] "L'Ecuyer-CMRG" "Inversion"     "Rejection"
```

### 2. Pre-computed RNG Streams per Task

Before any tasks execute, bayesim pre-computes an independent RNG stream
for each task:

``` r
# Create 5 independent RNG streams from seed 42
# Note: create_task_rng_streams() is an internal API used here for demonstration
streams <- bayesim:::create_task_rng_streams(42, 5)

# Each stream is a complete RNG state
str(streams[[1]])
#>  int [1:7] 10407 -2133391687 507561766 1260545903 1362917092 -1772566379 -1344458670

# These streams are independent---advancing one doesn't affect others
# This is important for reproducibility regardless of execution order
```

Each task receives its own stream via the `rng_seed` column in the task
grid. A scalar task seed is also derived from that stream and passed
into the data generator and fitter context. This means:

- Task execution order does not matter—each task uses its pre-assigned
  stream
- Parallel execution produces identical results to sequential execution
- Resuming from a checkpoint produces identical results to an
  uninterrupted run

### 3. Deterministic Task Ordering

Tasks are always processed in a deterministic, lexicographic order based
on their task IDs:

``` r
# Task IDs follow the pattern: d..._f..._r...
# - d...: zero-padded data grid index
# - f...: zero-padded fit grid index
# - r...: zero-padded replicate number
# Padding width is chosen from the study size so ordering stays stable.

# Example: With 2 data configs, 2 fit configs, and 3 replicates:
# d001_f001_r00001, d001_f001_r00002, d001_f001_r00003,
# d001_f002_r00001, d001_f002_r00002, d001_f002_r00003,
# d002_f001_r00001, d002_f001_r00002, d002_f001_r00003,
# d002_f002_r00001, d002_f002_r00002, d002_f002_r00003
```

This deterministic ordering ensures that:

- Tasks are always processed in the same sequence
- Checkpoint resume always continues from the same point
- Results are assembled in a predictable order

### 4. Configuration Fingerprinting

Every simulation configuration gets a unique cryptographic fingerprint:

``` r
# Create a simple configuration
config <- simulation_config(
  data_grid = data.frame(n = c(100, 200)),
  fit_grid = data.frame(model = "linear"),
  data_generator = function(data_spec, seed, task_ctx) {
    # Note: seed is a scalar task seed.
    # The simulation engine also restores the task RNG stream before each call.
    n <- data_spec$n
    x <- rnorm(n)
    y <- 1 + 2 * x + rnorm(n)
    list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(intercept = 1, slope = 2),
      vars_of_interest = c("intercept", "slope"),
      references = c(intercept = 0, slope = 0),
      meta = list()
    )
  },
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 3L,
  seed = 42L
)

# Compute the configuration fingerprint
fingerprint <- compute_config_fingerprint(config)
print(fingerprint)
#> [1] "33e998951dcaad6ce989f0669adf3820b382be6702c60d640c79869a59fa6187"
```

The fingerprint is a SHA256 hash of the normalized configuration,
excluding only runtime-specific fields (`result_path`,
`checkpoint_every`). This fingerprint:

- Uniquely identifies a simulation configuration
- Enables validation during checkpoint resume
- Allows detection of configuration changes

## Reproducibility Across Execution Modes

### Sequential vs Parallel Execution

Supported sequential and parallel executions use the same pre-assigned
task streams, so results stay aligned when the R, package, and
model-backend environment matches:

``` r
# Sequential execution
config <- simulation_config(..., seed = 42L)
result_seq <- run_simulation(config)

# Parallel execution
future::plan(future::multisession, workers = 4)
result_par <- run_simulation(config)

# These will produce identical results:
identical(result_seq$summary, result_par$summary)
```

This is possible because each task has its own pre-computed RNG stream.
Tasks can execute in any order without affecting their assigned
randomness.

### Resume from Checkpoint

Checkpoint/resume preserves the same task streams and completed-task
ledger:

``` r
# First run (interrupted)
config <- simulation_config(
  ...,
  seed = 42L,
  result_path = "my_simulation",
  checkpoint_every = 10L
)
result1 <- run_simulation(config)
# Suppose this gets interrupted after 50 tasks

# Resume from checkpoint
result2 <- run_simulation(config, resume = "must")

# result2 matches the uninterrupted run under the same software/backend stack
```

The resume mechanism:

1.  Loads the checkpoint ledger (task grid with statuses)
2.  Validates the configuration fingerprint matches
3.  Skips already-completed tasks
4.  Continues with pending tasks using their pre-assigned RNG streams
5.  Reassembles the same final results under the same software/backend
    stack

### Different Machines and Operating Systems

Simulations are reproducible across different machines and operating
systems **when**:

- The R version matches (see below for details)
- The bayesim version matches
- All package dependencies match
- The configuration and seed are identical

This allows you to:

- Start a simulation on your laptop and resume on a cluster
- Share results with collaborators who use different operating systems
- Archive simulations for long-term reproducibility

## What Is Deterministic

bayesim guarantees the following are reproducible within a matched
software/backend environment:

### Same Seed Produces Same Results

Given identical configuration and seed, you get identical results:

``` r
# Create a simple simulation configuration
data_gen <- function(data_spec, seed, task_ctx) {
  # Note: seed is a scalar task seed.
  # The simulation engine also restores the task RNG stream before each call.
  n <- data_spec$n
  x <- rnorm(n)
  y <- 2 + 3 * x + rnorm(n, sd = 0.5)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(intercept = 2, slope = 3, sigma = 0.5),
    vars_of_interest = c("intercept", "slope"),
    references = c(intercept = 0, slope = 0),
    meta = list()
  )
}

config <- simulation_config(
  data_grid = data.frame(n = 50, sigma = 0.5),
  fit_grid = data.frame(model = "test"),
  data_generator = data_gen,
  fitter = MockFitter(),
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 5L,
  seed = 123L
)

# Run twice with same seed
result_a <- run_simulation(config, progress = FALSE)
#> ℹ Starting simulation with 5 tasks
result_b <- run_simulation(config, progress = FALSE)
#> ℹ Starting simulation with 5 tasks

# Results are identical
all.equal(result_a$summary$rmse__value, result_b$summary$rmse__value)
#> [1] TRUE
all.equal(result_a$summary$bias__value, result_b$summary$bias__value)
#> [1] TRUE
```

### Resume Produces Same Final Results

Resuming from a checkpoint produces the same final results as an
uninterrupted run:

``` r
# This is demonstrated conceptually - in practice you would need
# to interrupt and resume an actual run

# Uninterrupted run
config <- simulation_config(
  ...,
  seed = 42L,
  result_path = "uninterrupted_run",
  checkpoint_every = 1000L  # Won't trigger during run
)
result_full <- run_simulation(config)

# Interrupted and resumed run
config2 <- simulation_config(
  ...,
  seed = 42L,
  result_path = "resumed_run",
  checkpoint_every = 10L
)
result_resumed <- run_simulation(config2, resume = "must")

# Final results are identical
identical(result_full$summary, result_resumed$summary)
```

### Task Execution Order

Tasks are planned in a deterministic order:

1.  Sorted lexicographically by task ID
2.  Each task uses its pre-assigned RNG stream
3.  Results are assembled in task ID order

## What is NOT Guaranteed

While bayesim strives for complete reproducibility, some aspects are
inherently non-deterministic:

### Exact Timing Values

Task execution times will vary between runs:

``` r
# Run the same simulation twice
config <- simulation_config(
  data_grid = data.frame(n = 100),
  fit_grid = data.frame(model = "test"),
  data_generator = data_gen,
  fitter = MockFitter(),
  metrics = list(rmse_metric()),
  n_replicates = 3L,
  seed = 42L
)

result1 <- run_simulation(config, progress = FALSE)
#> ℹ Starting simulation with 3 tasks
result2 <- run_simulation(config, progress = FALSE)
#> ℹ Starting simulation with 3 tasks

# Timing will differ
timing1 <- result1$timing$total
timing2 <- result2$timing$total
cat("Run 1:", timing1, "seconds\n")
#> Run 1: 0.04330945 seconds
cat("Run 2:", timing2, "seconds\n")
#> Run 2: 0.04285169 seconds
```

This is expected and doesn’t affect the scientific validity of results.

### Order of Warning Messages

If multiple tasks produce warnings, the order in which they appear may
vary:

- In sequential execution, warnings appear in task order
- In parallel execution, warning order depends on task completion order
- The warnings themselves are deterministic, just their order may vary

### Results Across R Versions

**Critical**: Results may differ across different versions of R, even
with the same seed. This is because:

- R’s random number generation may change between versions
- Underlying algorithms (e.g., for linear algebra) may be updated
- Floating-point behavior can vary

**Best practice**: Always record the R version used for your
simulations.

### Results Across Package Versions

Similarly, results may differ across versions of bayesim or its
dependencies:

- Algorithm improvements may change results
- Bug fixes may alter behavior
- New features may change defaults

**Best practice**: Use package version pinning (e.g., renv, packrat) for
critical simulations.

## Best Practices

### Always Set a Seed

Always explicitly set a seed in your simulation configuration:

``` r
# Good: explicit seed
config <- simulation_config(
  ...,
  seed = 42L  # Use any integer
)

# Without a seed, your results won't be reproducible
```

Choose a seed arbitrarily (many researchers use their birthday, favorite
number, or random.org). The specific value doesn’t matter—only that it’s
documented.

### Version Control Your Configuration

Store your simulation configuration in version control:

``` r
# In your R script or R Markdown document
simulation_config(
  data_grid = read.csv("data_grid.csv"),
  fit_grid = read.csv("fit_grid.csv"),
  data_generator = my_data_gen,  # Defined in the same file
  fitter = my_fitter,            # Defined in the same file
  metrics = list(rmse_metric(), bias_metric()),
  n_replicates = 1000L,
  seed = 42L
)
```

This ensures that:

- The exact configuration is preserved
- Changes are tracked over time
- Others can reproduce your setup

### Use Checkpointing for Long Runs

For long-running simulations, always enable checkpointing:

``` r
config <- simulation_config(
  ...,
  result_path = "path/to/results",  # Enable checkpointing
  checkpoint_every = 50L             # Save every 50 tasks
)

# Run with resume capability
result <- run_simulation(config, resume = "auto")
```

This protects against:

- System crashes
- Power failures
- Time limits on clusters
- Accidental interruption

### Record R and Package Versions

Document your computational environment:

``` r
# Record R version
R.version.string
#> [1] "R version 4.5.3 (2026-03-11)"

# Record package versions
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> Random number generation:
#>  RNG:     L'Ecuyer-CMRG 
#>  Normal:  Inversion 
#>  Sample:  Rejection 
#>  
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] future_1.70.0 bayesim_1.0.1
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.1         cli_3.6.5           knitr_1.51         
#>  [4] rlang_1.1.7         xfun_0.56           generics_0.1.4     
#>  [7] textshaping_1.0.5   S7_0.2.1            jsonlite_2.0.0     
#> [10] glue_1.8.0          listenv_0.10.1      future.apply_1.20.2
#> [13] htmltools_0.5.9     ragg_1.5.1          sass_0.4.10        
#> [16] rmarkdown_2.30      evaluate_1.0.5      jquerylib_0.1.4    
#> [19] tibble_3.3.1        fastmap_1.2.0       yaml_2.3.12        
#> [22] lifecycle_1.0.5     compiler_4.5.3      dplyr_1.2.0        
#> [25] codetools_0.2-20    fs_1.6.7            pkgconfig_2.0.3    
#> [28] systemfonts_1.3.2   digest_0.6.39       R6_2.6.1           
#> [31] tidyselect_1.2.1    parallelly_1.46.1   pillar_1.11.1      
#> [34] parallel_4.5.3      magrittr_2.0.4      bslib_0.10.0       
#> [37] tools_4.5.3         withr_3.0.2         globals_0.19.1     
#> [40] pkgdown_2.2.0       cachem_1.1.0        desc_1.4.3
```

For publication, include a section like:

> “Simulations were conducted using R version 4.3.1 with bayesim version
> 0.1.0. The computational environment is available at
> \[repository/link\].”

### Environment Management

For stronger reproducibility guarantees, use environment management
tools to ensure consistent package versions and R environments across
machines and sessions.

#### renv: Project-Local Dependencies

**renv** creates a project-specific package library with a lockfile
capturing exact package versions:

``` r
# Initialize renv in your project
renv::init()

# Snapshot current package versions to lockfile
renv::snapshot()

# Restore packages from lockfile (on another machine or session)
renv::restore()
```

The `renv.lock` file ensures collaborators and future-you install
identical package versions. This is necessary for long-term
reproducibility of simulation results.

#### rv and Containers: Full Environment Capture

For capturing the complete R environment (R version + packages + system
libraries):

- **rv** is a newer tool for managing R environments with
  version-controlled snapshots
- **Docker/containers** provide the strongest guarantees by capturing
  the entire OS, R version, and all dependencies

``` r
# Example: Using rv to manage environments
rv::init()
rv::snapshot()
rv::restore()
```

Container-based approaches (Docker, Singularity) are recommended for
published research where reproducibility must survive years of software
evolution.

## Verifying Reproducibility

Here’s a complete example demonstrating reproducibility:

``` r
library(bayesim)

# Define a data generator
my_data_generator <- function(data_spec, seed, task_ctx) {
  # Note: seed is a scalar task seed.
  # The simulation engine also restores the task RNG stream before each call.
  n <- data_spec$n
  
  # Generate synthetic data with known parameters
  x <- rnorm(n)
  y <- data_spec$intercept + data_spec$slope * x + rnorm(n, sd = data_spec$sigma)
  
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(
      intercept = data_spec$intercept,
      slope = data_spec$slope,
      sigma = data_spec$sigma
    ),
    vars_of_interest = c("intercept", "slope"),
    references = c(intercept = 0, slope = 0),
    meta = list()
  )
}

# Create configuration
config <- simulation_config(
  data_grid = data.frame(
    n = c(50, 100),
    intercept = 1,
    slope = 2,
    sigma = 0.5
  ),
  fit_grid = data.frame(
    model = "linear"
  ),
  data_generator = my_data_generator,
  fitter = MockFitter(),
  metrics = list(
    rmse_metric(),
    bias_metric()
  ),
  n_replicates = 10L,
  seed = 42L
)

# Run simulation twice with same seed
cat("Running simulation for the first time...\n")
#> Running simulation for the first time...
result1 <- run_simulation(config, progress = FALSE)
#> ℹ Starting simulation with 20 tasks

cat("Running simulation for the second time...\n")
#> Running simulation for the second time...
result2 <- run_simulation(config, progress = FALSE)
#> ℹ Starting simulation with 20 tasks

# Verify results are identical
cat("\n=== Reproducibility Check ===\n")
#> 
#> === Reproducibility Check ===

# Check summary structure
structures_match <- identical(names(result1$summary), names(result2$summary))
cat("Summary structures match:", structures_match, "\n")
#> Summary structures match: TRUE

# Check metric values
rmse_match <- all.equal(result1$summary$rmse__value, result2$summary$rmse__value)
bias_match <- all.equal(result1$summary$bias__value, result2$summary$bias__value)
cat("RMSE values match:", isTRUE(rmse_match), "\n")
#> RMSE values match: TRUE
cat("Bias values match:", isTRUE(bias_match), "\n")
#> Bias values match: TRUE

# Check task statuses
status_match <- identical(result1$task_grid$status, result2$task_grid$status)
cat("Task statuses match:", status_match, "\n")
#> Task statuses match: TRUE

# Overall reproducibility
is_reproducible <- structures_match && isTRUE(rmse_match) && 
                   isTRUE(bias_match) && status_match
cat("\nOverall reproducibility:", is_reproducible, "\n")
#> 
#> Overall reproducibility: TRUE
```

### Testing with Different Seeds

To verify that different seeds produce different (but still
reproducible) results:

``` r
# Same configuration, different seeds
config1 <- simulation_config(
  data_grid = data.frame(n = 100, intercept = 1, slope = 2, sigma = 0.5),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_data_generator,
  fitter = MockFitter(),
  metrics = list(rmse_metric()),
  n_replicates = 5L,
  seed = 42L
)

config2 <- simulation_config(
  data_grid = data.frame(n = 100, intercept = 1, slope = 2, sigma = 0.5),
  fit_grid = data.frame(model = "linear"),
  data_generator = my_data_generator,
  fitter = MockFitter(),
  metrics = list(rmse_metric()),
  n_replicates = 5L,
  seed = 123L  # Different seed
)

result1 <- run_simulation(config1, progress = FALSE)
#> ℹ Starting simulation with 5 tasks
result2 <- run_simulation(config2, progress = FALSE)
#> ℹ Starting simulation with 5 tasks

# Results should be different (different random data)
cat("Same seed produces same results: ",
    identical(result1$summary$rmse__value, result2$summary$rmse__value), "\n")
#> Same seed produces same results:  TRUE
cat("Different seeds produce different results: ",
    !identical(result1$summary$rmse__value, result2$summary$rmse__value), "\n")
#> Different seeds produce different results:  FALSE
```

## Summary

bayesim provides strong determinism guarantees for simulation studies:

| Aspect                          | Guarantee                                 |
|---------------------------------|-------------------------------------------|
| Same seed → same results        | ✅ Guaranteed                             |
| Sequential vs parallel          | ✅ Identical in supported execution modes |
| Resume from checkpoint          | ✅ Identical to uninterrupted run         |
| Task execution order            | ✅ Deterministic                          |
| Cross-platform (same R version) | ⚠️ Expected under controlled environment  |
| Timing values                   | ❌ Not guaranteed                         |
| Warning message order           | ❌ Not guaranteed                         |
| Across R versions               | ❌ Not guaranteed                         |
| Across package versions         | ❌ Not guaranteed                         |

By following the best practices outlined in this vignette—always setting
a seed, version controlling your configuration, using checkpointing,
recording your computational environment, and using environment
management tools—you ensure that your simulation studies remain
reproducible and scientifically valid.

For more information, see:

- [`vignette("getting-started")`](https://sims1253.github.io/bayesim/articles/getting-started.md)
  for basic bayesim usage
- [`vignette("simulation-study")`](https://sims1253.github.io/bayesim/articles/simulation-study.md)
  for a complete example
- [`vignette("memory-management")`](https://sims1253.github.io/bayesim/articles/memory-management.md)
  for handling large simulations
