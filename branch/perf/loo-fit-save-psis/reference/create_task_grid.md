# Create Task Grid from Configuration

Generates the complete, deterministic task table for a simulation. Tasks
are ordered lexicographically by task_id for reproducibility.

## Usage

``` r
create_task_grid(config)
```

## Arguments

- config:

  A SimulationConfig S7 object.

## Value

A tibble with one row per task, containing columns:

- `task_id`: Character identifier in format "dXXX_fXXX_rXXXXX"

- `data_idx`: Integer index into the data_grid (1-based)

- `fit_idx`: Integer index into the fit_grid (1-based)

- `rep_idx`: Integer replicate number (1 to n_replicates)

- `rng_seed`: List column containing precomputed RNG stream for each
  task

- `status`: Character status, initialized to "pending"

## Examples

``` r
if (FALSE) { # \dontrun{
config <- simulation_config(
  data_grid = data.frame(n = c(100, 500)),
  fit_grid = data.frame(model = c("baseline", "full")),
  data_generator = my_data_gen,
  n_replicates = 10L,
  seed = 42L
)

task_grid <- create_task_grid(config)
# task_grid has 2 * 2 * 10 = 40 rows
} # }
```
