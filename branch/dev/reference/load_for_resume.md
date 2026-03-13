# Load Run for Resume

Loads previous run state for resumption. Validates schema compatibility
and configuration fingerprint before allowing resume.

## Usage

``` r
load_for_resume(result_path, config)
```

## Arguments

- result_path:

  Character; path to results directory containing checkpoints.

- config:

  SimulationConfig; current configuration object.

## Value

A list with elements:

- `task_grid`: Task grid with restored status from checkpoint

- `prior_results`: Data frame of results from checkpoint

- `checkpoint_id`: ID of the checkpoint being resumed from

## Details

The function performs the following validation steps:

1.  Reads and validates run_manifest.json

2.  Checks schema version compatibility

3.  Computes and compares configuration fingerprints

4.  Finds the most recent valid checkpoint

5.  Rebuilds task grid with status from checkpoint

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(...)
resume_state <- load_for_resume("/path/to/results", config)
task_grid <- resume_state$task_grid
prior_results <- resume_state$prior_results
} # }
```
