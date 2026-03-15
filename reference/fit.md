# Fit a Bayesian Model

Main fitting method for Bayesian model fitters.

## Usage

``` r
fit(fitter, data_bundle, fit_spec, seed, task_ctx)
```

## Arguments

- fitter:

  An S7 Fitter object

- data_bundle:

  A list containing:

  - `train`: Training data (data.frame)

  - `test`: Test data (data.frame, may be NULL)

  - `response`: Name of response variable (character)

  - `true_params`: True parameter values (named numeric, may be NULL)

  - Additional model-specific data

- fit_spec:

  A list or single-row data.frame from fit_grid containing:

  - Model specification parameters

  - Prior specifications

  - Algorithm settings

- seed:

  Integer seed for reproducibility

- task_ctx:

  A list with simulation context:

  - `task_id`: Unique task identifier

  - `data_idx`: Index of data generation settings

  - `fit_idx`: Index of fitting settings

  - `rep_idx`: Replication index

## Value

A `bayesim_fit_result` S3 object containing:

- success: TRUE/FALSE

- fit: backend object

- draws: matrix with colnames or NULL

- diagnostics: named list

- timing: list with total, warmup, sample

- warnings: character vector

- error: NULL or condition
