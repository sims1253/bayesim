# bayesim

bayesim runs simulation studies of Bayesian models in R. Define how to
generate data, which models to fit, and what to measure. It runs the
replicates and summarizes estimation error, interval coverage, or
calibration. Parallel execution and checkpoints support longer studies.

## Install

Install the development version from GitHub:

``` r

# install.packages("pak")
pak::pak("sims1253/bayesim")
```

The example below runs without Stan. For brms and CmdStan setup, see
[brms
studies](https://sims1253.github.io/bayesim/articles/brms-studies.html)
and
[fitters](https://sims1253.github.io/bayesim/articles/custom-fitters.html).

## Try a study

Does a 90% credible interval cover the true slope about 90% of the time?
Compare three sample sizes using exact Bayesian linear regression:

``` r

library(bayesim)

gen <- function(data_spec, task_ctx) {
  x <- rnorm(data_spec$n)
  y <- 1 + 0.5 * x + rnorm(data_spec$n)
  list(
    train = data.frame(y = y, x = x),
    response = "y",
    true_params = c(x = 0.5),
    vars_of_interest = "x"
  )
}

config <- simulation_config(
  data_grid = data.frame(n = c(20, 50, 200)),
  fit_grid = data.frame(model = "linear"),
  data_generator = gen,
  fitter = LinearRegressionFitter(),
  metrics = list(coverage_metric(prob = 0.90)),
  n_replicates = 100L,
  seed = 42L
)

result <- run_simulation(config, progress = FALSE, verbose = FALSE)
summary <- summarize_simulation(result, by = "data_n", metrics = "coverage__by_param__x")
```

| Sample size | Coverage |  MCSE |
|------------:|---------:|------:|
|          20 |     0.77 | 0.042 |
|          50 |     0.88 | 0.032 |
|         200 |     0.95 | 0.022 |

This runs 300 fits: 3 sample sizes × 1 model × 100 replicates. Compare
the coverage with 0.90; its Monte Carlo standard error (MCSE) describes
uncertainty from the finite number of replicates. Use more replicates
for a precise estimate. Add `workers = 2` to
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md)
to run in parallel.

Here the smallest sample has coverage below 0.90. Exact posterior
inference does not guarantee nominal interval coverage at a fixed
parameter value; this is what the study measures.

## Guides

- [Getting
  started](https://sims1253.github.io/bayesim/articles/getting-started.html):
  build a study and interpret its results.
- [Study
  design](https://sims1253.github.io/bayesim/articles/design-of-simulation-studies.html)
  and
  [SBC](https://sims1253.github.io/bayesim/articles/sbc-and-calibration.html):
  choose what to measure.
- [Parallel runs and
  checkpoints](https://sims1253.github.io/bayesim/articles/parallel-and-hpc.html):
  run larger studies and resume interrupted work.
- [Reproducibility](https://sims1253.github.io/bayesim/articles/reproducibility.html):
  seeds, fingerprints, and their limits.
- [Function
  reference](https://sims1253.github.io/bayesim/reference/index.html):
  fitters, generators, metrics, and analysis tools.

Report bugs or request features in the [issue
tracker](https://github.com/sims1253/bayesim/issues).
