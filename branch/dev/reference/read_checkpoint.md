# Read checkpoint

Reads a checkpoint by ID, or the latest valid checkpoint if ID is not
specified. Verifies checksums before returning data.

## Usage

``` r
read_checkpoint(result_path, checkpoint_id = NULL)
```

## Arguments

- result_path:

  Character string giving the base path for results. If NULL, the
  function returns NULL immediately.

- checkpoint_id:

  Integer specifying the checkpoint ID to read. If NULL (default), reads
  the latest checkpoint from latest.json.

## Value

A list with the following elements, or NULL if checkpoint not found or
invalid:

- `checkpoint_id` - Integer ID of the checkpoint

- `meta` - List of checkpoint metadata

- `task_grid` - Tibble/data.frame of task grid with status

- `results_df` - Data frame of task results

## Details

The function performs the following steps:

1.  If checkpoint_id is NULL, read latest.json to get the latest ID

2.  Construct checkpoint directory path from the ID

3.  Verify directory exists

4.  Verify checksums match

5.  Read meta.json, ledger.rds, and results.rds

6.  Return assembled checkpoint data

If the checkpoint has invalid checksums, a warning is issued and NULL is
returned. This allows the caller to fall back to earlier checkpoints.

## See also

[`write_checkpoint()`](https://sims1253.github.io/bayesim/reference/write_checkpoint.md),
[`list_checkpoints()`](https://sims1253.github.io/bayesim/reference/list_checkpoints.md),
[`validate_checkpoint_fingerprint()`](https://sims1253.github.io/bayesim/reference/validate_checkpoint_fingerprint.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Read latest checkpoint
checkpoint <- read_checkpoint("/path/to/results")

# Read specific checkpoint
checkpoint <- read_checkpoint("/path/to/results", checkpoint_id = 5)
} # }
```
