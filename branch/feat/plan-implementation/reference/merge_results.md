# Merge Prior and New Results

Combines results from a resumed run with new results from continued
execution. Duplicate task rows must be identical or an integrity error
is raised.

## Usage

``` r
merge_results(prior_results, new_results)
```

## Arguments

- prior_results:

  Data frame of results from checkpoint.

- new_results:

  Data frame of results from new execution.

## Value

Combined data frame with unique task_id values.

## Details

If both inputs are NULL or empty, returns NULL. If only one input has
data, returns that input. Otherwise, removes any task_ids from
prior_results that appear in new_results, then combines.
