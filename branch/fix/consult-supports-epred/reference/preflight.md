# Preflight check for a simulation configuration

Inspects a configuration *before* a long run and reports: total task
count, grid dimensions, metrics and their `needs` vs the fitter's
capabilities (surfacing mismatch warnings before a 10-hour run), whether
mirai daemons are set, and (for BrmsFitter) the number of distinct
models to compile.

When run automatically at the top of
[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md),
it prints a condensed one-line summary.

## Usage

``` r
preflight(config, pilot = FALSE, condensed = FALSE)
```

## Arguments

- config:

  A `SimulationConfig`.

- pilot:

  Logical; if TRUE, run a single pilot task and extrapolate its
  wall-clock time to an estimate of the total study time (default
  FALSE). The estimate is printed and returned in the `pilot_seconds`
  and `estimated_total_seconds` elements of the returned list. A pilot
  that cannot execute (e.g. a Stan fitter with no CmdStan installed)
  yields NA estimates rather than failing the preflight. Note that for
  Stan fitters a pilot compiles the model, so it is not instant.

- condensed:

  Logical; if TRUE, print a single one-line summary instead of the full
  report (used by run_simulation).

## Value

Invisibly, a named list of preflight information.

## See also

[`run_simulation()`](https://sims1253.github.io/bayesim/reference/run_simulation.md),
[`simulation_config()`](https://sims1253.github.io/bayesim/reference/simulation_config.md)

## Examples

``` r
if (FALSE) { # \dontrun{
preflight(config)
} # }
```
