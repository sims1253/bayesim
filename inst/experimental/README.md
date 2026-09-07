# Experimental study grammar

This directory contains an executable interface experiment. Load it with:

```r
source(system.file("experimental", "study-grammar.R", package = "bayesim"))
```

It is separate from the supported `simulation_config()` engine. Its purpose is
to test whether real studies can express their design, artifact needs and
comparison rules without writing an execution loop.

The case-study vignettes build plans; they do not run the original simulations.
The small analytic example exercises the reference runner.

## Function contracts

- `study(name, conditions, generate, prepare = NULL, version = "1")` starts a
  declaration. `conditions` has a unique character `condition_id` column.
  `generate(condition, context)` returns any R object. Optional
  `prepare(condition, context)` runs once per condition; its output is available
  as `context$prepared`. Increment `version` when external generation or
  preparation dependencies change.
- `study_method(fit, extract = list(), settings = list(), version = "1")`
  describes a method. `fit(data, context)` returns its native result. Each named
  extractor is `function(fit, data, context)` and returns an artifact, such as
  chain-preserving draws, diagnostics or pointwise predictive scores. Method
  settings are in `context$settings`. Use namespace-qualified calls on workers.
- `with_methods(study, ...)` adds named methods. Names are persistent IDs.
  `do.call(with_methods, c(list(study), named_methods))` adds a programmatically
  constructed collection.
- `with_measure(study, name, compute, needs = character(), version = "1")`
  declares a per-fit measurement. `compute(artifacts, context)` returns a data
  frame, usually with `target` and `value`. List columns are allowed. `needs`
  names extractors, or the built-in artifacts `data` and `fit`.
- `with_comparison(study, name, compute, by = "dataset_id", needs = character(),
  version = "1")` declares a grouped calculation. `compute(results, artifacts,
  context)` gets the group's measurement rows and a list of artifact bundles
  named by method ID. Artifact-consuming groups must stay within one dataset;
  wider groups can consume saved measurement rows. All method attempts remain
  in `context$attempts`, including failures. Return a data frame.
- `with_retention(study, artifacts = character(), compress = "gzip")` selects
  artifact names to save. Measurements, errors and identity metadata are always
  saved. Retaining `fit` is explicit. Extractors select parameters and preserve
  chain information; retention does not convert draw types.
- `plan_study(study, replicates, seed = 1L)` validates declarations and describes
  counts, required artifacts, computation stages and retention. It does not
  generate data, prepare models, inspect an existing cache or estimate sizes
  without data.
- `run_study(study, replicates, seed = 1L, path = NULL, workers = 1L)` executes
  the plan. `evaluate_replicate(study, condition_id, replicate, seed = 1L,
  path = NULL)` is the same dataset-level evaluator for external schedulers.
- `assess_study(run, policy = NULL)` applies an optional
  `function(measurements, attempts)` policy to select measurement rows, then
  reruns comparisons that depend only on those rows. It retains raw rows and
  an exclusion table. Artifact-consuming comparisons are computed during
  execution and are not silently reinterpreted under a later selection policy.

`context` contains study name, condition, dataset ID, replicate, method ID,
settings, prepared output and a stage-specific integer seed. The ambient RNG
also follows that seed. Measurement and comparison streams are independent of
fitting and generation streams.
