# Experimental case studies

The study grammar that started here is now part of the package: `study()`,
`study_method()`, `with_methods()`, `with_measure()`, `with_comparison()`,
`with_retention()`, `plan_study()`, `run_study()`, `evaluate_replicate()`,
`assess_study()` and `summarise_study()`. Load it with `library(bayesim)`.
See `?study` and `vignette("experimental-study-grammar", package = "bayesim")`.

This directory keeps one file:

- `likelihood-case.R` declares the likelihood comparison study of
  `vignette("experimental-likelihood-study", package = "bayesim")`. Load it
  after `library(bayesim)` with

  ```r
  source(system.file("experimental", "likelihood-case.R", package = "bayesim"))
  ```

## Changes from the experimental runner

- The store is content-addressed: one file per fit, named by a digest of
  everything that determines it. There is no manifest. A store can hold
  several seeds and several versions of a study.
- Several R processes can run the same study on one store at the same time.
- `run_study()` takes `backend = mirai_backend(n)` instead of `workers = n`.
  Jobs start as soon as a worker is free; there are no batches.
- Store keys include the values and helper functions that a function uses,
  not only its body and version. `plan_study()` stops when a function uses a
  value that cannot be sent to a worker, or one that is too large.
- `prepare()` runs as a job, on a worker, and its value is stored.
- A failed fit is fitted again with `run_study(retry = TRUE)`.
- A condition field or method setting named `id` is an error (issue #110).
