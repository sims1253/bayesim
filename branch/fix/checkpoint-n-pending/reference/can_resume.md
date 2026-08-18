# Check if Resumable Run Exists

Checks whether a valid checkpoint exists that can be resumed. A run can
be resumed if:

- Both run_manifest.json and latest.json exist

- latest.json points to a valid checkpoint_id

- The referenced checkpoint can be read and validated

## Usage

``` r
can_resume(result_path)
```

## Arguments

- result_path:

  Character; path to results directory containing checkpoints.

## Value

TRUE if a valid run can be resumed, FALSE otherwise.

## Details

This is a cheap existence/validity probe: the checkpoint is validated
(checksums, ledger, shard integrity) with `load_outcomes = FALSE`, so
the full outcome history is never deserialized here. Callers that need
the outcomes use
[`load_for_resume()`](https://sims1253.github.io/bayesim/reference/load_for_resume.md)
instead.
