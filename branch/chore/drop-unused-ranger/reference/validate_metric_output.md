# Validate Metric Output

Validates that metric output conforms to the expected schema. Metric
output must be a named list where all names are non-empty strings and
values are either scalar atomic types or named numeric vectors.

## Usage

``` r
validate_metric_output(output, metric_name)
```

## Arguments

- output:

  The output list to validate.

- metric_name:

  Character string identifying the metric (for error messages).

## Value

Invisible `output` if validation passes. Otherwise, an error is raised.

## Validation Rules

- output must be a list

- all elements must have non-empty names

- values must be one of:

  - scalar logical, integer, double, or character

  - named numeric vector

- no NULL values allowed

- no nested lists, data frames, or matrices allowed

## Examples
