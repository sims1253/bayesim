# Bayesim

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)
[![R-CMD-check](https://github.com/sims1253/bayesim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/R-CMD-check.yaml)
[![Tests](https://github.com/sims1253/bayesim/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/test-coverage.yaml)
[![Codecov test
coverage](https://codecov.io/gh/sims1253/bayesim/graph/badge.svg)](https://app.codecov.io/gh/sims1253/bayesim)
[![GH-Pages](https://github.com/sims1253/bayesim/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/sims1253/bayesim/actions/workflows/pkgdown.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

A simulation framework for Bayesian models based on [brms](https://github.com/paul-buerkner/brms/).

The main function is `full_simulation` with the main arguments being `data_gen_confs`, `data_gen_fun`, `fit_confs` and `metrics`.
Bayesim will generate datasets by passing `data_gen_confs` rows to `data_gen_fun` and fit each model defined by `fit_confs` on each generated dataset. It then calculates all of the defined `metrics` for each model. This is, as of now, done in a fully crossed fashion.

## Define Data Simulation

Data simulation consists of two parts. A `data_gen_fun` function and a `data_gen_confs` dataframe. Bayesim will feed each row of `data_gen_confs` into `data_gen_fun` to generate each individual dataset.

The only strictly necessary columns in `data_gen_confs` are:

- `dataset_N`, the number of datasets that should be simulated per configuration. This is what Bayesim parallelize over.

- `id`, a unique identifier string that is used to save the results per configuration row.

- `vars_of_interest`, if you want metrics calculated for individual parameters. It should be a list of variable names used in `data_gen_fun`.

`data_gen_fun` should output a named list that contains the following parts:

- `dataset`, a dataframe that is fed into `brms::brm`

- `testing_data`, a dataframe that is used as `newdata` argument for certain out-of-sample metrics.

- `data_gen_output`, a named list that contains all other information that you want to export from the data generating function. This should usually include all the input arguments (see example for how to get those easily), and a `references` list, that contains the reference values for all `vars_of_interest` variables (again, see the example for how to easily get those).
 An example is presented below:

```{r}
constant_linpred_dgp <- function(data_N,
                                 data_link,
                                 data_family,
                                 seed = NULL,
                                 testing_data = TRUE,
                                 vars_of_interest = list("mu"),
                                 mean = 0,
                                 ...) {
  arguments <- as.list(c(as.list(environment()), list(...)))
  arguments$seed <- NULL

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (testing_data) {
    data_gen_size <- data_N * 2
  } else {
    data_gen_size <- data_N
  }
  dataset <- data.frame()
  mu = rnorm(n = 1, mean = x, sd = 1)
  y = rnorm(n = data_gen_size, mean = mu, sd = 1)

  # This creates a list of values for each of the vars_of_interest. 
  arguments$references <- lapply(
    unlist(vars_of_interest),
    function(x) get(x)
  )

  data_gen_output <- list()
   # Anything in addition to the function arguments you want to save about
   # the data generation process ie. If you are resampling the number
   # of invalid samples.
  )
  data_gen_output <- c(data_gen_output, arguments)

  if (testing_data) {
    return(
      list(
        dataset = list(y = dataset[1:data_N, ]),
        testing_data = list(y = dataset[(data_N + 1):data_gen_size, ]),
        data_gen_output = data_gen_output
      )
    )
  } else {
    return(
      list(
        dataset = dataset,
        testing_data = NULL,
        data_gen_output = data_gen_output
      )
    )
  }
}

```

## Define Fit Configurations

Fit configurations currently are dataframes with the following columns:

- `fit_family`, see [brms_family_lookup](R/ll_lookup.R#L10) for supported families.

- `fit_link`, see [link_lookup](R/ll_lookup.R#L70) for supported families.

- `formula`, a string that allows conversion via `brms::brmsformula`

- `prior`, gets passed to `brms::brm` directly.

## Define Metrics

Metrics are defined via a list of string identifiers. The supported metrics are:

### Variable summaries

`"v_mean"`
`"v_sd"`
`"v_median"`
`"v_mad"`
`"v_pos_prob"`
`"v_quantiles"`
`"v_bias"`
`"v_rmse"`
`"v_mae"`
`"v_mse"`
`"v_true_percentile"`

### Global MCMC Diagnostics

`"divergent_transitions_rel"`
`"divergent_transitions_abs"`
`"rstar"`
`"bad_pareto_ks"`
`"pareto_k_values"`
`"time_per_sample"`

### Variable MCMC Diagnostics

`"rhat"`
`"ess_bulk"`
`"ess_tail"`

### Predictive Metrics

`"elpd_loo"`
`"elpd_loo_pointwise"`
`"elpd_loo_pointwise_summary"`
`"elpd_test"`
`"elpd_test_pointwise_summary"`
`"rmse_loo"`
`"rmse_loo_pointwise"`
`"rmse_loo_pointwise_summary"`
`"rmse_test"`
`"rmse_test_pointwise_summary"`
`"r2_loo"`
`"r2_loo_pointwise"`
`"r2_loo_pointwise_summary"`
`"r2_test"`
`"r2_test_pointwise_summary"`

### Posterior sample based metrics

`"log_lik_pointwise"`
`"log_lik_summary"`
`"ppred_summary_y_scaled"`
`"ppred_pointwise"`
`"residuals"`
`"posterior_linpred"`
`"posterior_linpred_transformed"`

### Observations

`"y_pointwise"`
`"y_pointwise_z_scaled"`
`"y_summaries"`

### Data

`"data_gen"`

### Fits

`"fit_gen"`

Or see [metric_lookup](R/metric_lookup.R#L11) for all currently implemented metrics.

## Additional Arguments

`seed`, sets a seed that will result in the rest of the simulation happening deterministically, conditional on the seed. Allows for reproduction of individual results or the entire simulation run later on.

### Output

- `result_path = "./"`, The path where the result .RDS files should be saved.

- `debug = FALSE`, `TRUE` will save all intermediate results as .RDS files in the `result_path` directory to support debugging.

### Stan Options

`stan_pars` should be a named list that contains the following arguments:

- `warmup`, directly passed to `brms::brm`

- `iter`, directly passed to `brms::brm`

- `chains`, directly passed to `brms::brm`

- `init`, directly passed to `brms::brm`

- `backend = "rstan"`, directly passed to `brms::brm` We recommend rstan due to instabilities of cmdstanr on clusters.

- `cmdstan_path`, useful when working with cmdstan on a computing cluster where cmdstan might not be installed in the default location. Use the path to the main directory, eg `"~/.cmdstan/cmdstan-2.29.2"`.

- `cmdstan_write_path`, directory for cmdstan to write compiled model files to. This should not be a temporary directory as those might get cleaned up during the simulation run.

### Using multiple Cores

- `ncores_simulation = 1`, If set to more than `1`, Bayesim will parallelize across datasets within each row of `data_gen_confs` using the specified number of processes.

- `cluster_type = "PSOCK"`, Defines the type of cluster used by the `parallel` package. Windows requires `PSOCK` however 'FORK` can save quite some time due to the repeated cluster setup times.

## Related Work

Bayesim has been used in the following projects:

- [Prediction can be safely used as a proxy for explanation in causally consistent Bayesian generalized linear models](https://arxiv.org/abs/2210.06927) [![DOI](https://zenodo.org/badge/453991253.svg)](https://zenodo.org/badge/latestdoi/453991253)
