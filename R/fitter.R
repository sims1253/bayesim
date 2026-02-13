#' @title Abstract Fitter Class
#' @description
#' Abstract base class for Bayesian model fitters in bayesim.
#'
#' The Fitter class defines the interface that all model fitters must implement.
#' It provides a consistent API for fitting Bayesian models, extracting posterior
#' draws, generating predictions, computing log-likelihoods, and performing
#' model diagnostics.
#'
#' @section Properties:
#' \itemize{
#'   \item `name`: Character string identifying the fitter (e.g., "stan", "brms")
#'   \item `supports_predictions`: Logical indicating if the fitter supports predictions
#'   \item `supports_log_lik`: Logical indicating if the fitter supports log-likelihood computation
#'   \item `supports_loo`: Logical indicating if the fitter supports LOO-CV
#' }
#'
#' @section Methods:
#' The following S7 generics must be implemented by subclasses:
#' \describe{
#'   \item{`fit(fitter, data_bundle, fit_spec, seed, task_ctx)`}{Main fitting method}
#'   \item{`extract_draws(fitter, fit_result, variables = NULL)`}{Extract posterior draws}
#'   \item{`predict(fitter, fit_result, newdata = NULL, seed = NULL)`}{Generate predictions}
#'   \item{`log_lik(fitter, fit_result, newdata = NULL)`}{Pointwise log-likelihood}
#'   \item{`loo(fitter, fit_result)`}{LOO-CV computation}
#'   \item{`diagnostics(fitter, fit_result)`}{Extract fit diagnostics}
#' }
#'
#' @section Creating Custom Fitters:
#' To create a custom fitter, extend this class and implement methods for the
#' S7 generics: `fit()`, `extract_draws()`, `predict()`, `log_lik()`, `loo()`, `diagnostics()`.
#'
#' @return An S7 class object representing the abstract Fitter
#' @export
#' @seealso [MockFitter] for a simple example implementation
Fitter <- S7::new_class(
  "Fitter",
  abstract = TRUE,
  properties = list(
    name = S7::new_property(S7::class_character),
    supports_predictions = S7::new_property(S7::class_logical, default = TRUE),
    supports_log_lik = S7::new_property(S7::class_logical, default = TRUE),
    supports_loo = S7::new_property(S7::class_logical, default = TRUE)
  )
)

# =============================================================================
# S7 Generics for Fitter methods
# =============================================================================

#' @title Fit a Bayesian Model
#' @description
#' Main fitting method for Bayesian model fitters.
#'
#' @param fitter An S7 Fitter object
#' @param data_bundle A list containing:
#'   \itemize{
#'     \item `train`: Training data (data.frame)
#'     \item `test`: Test data (data.frame, may be NULL)
#'     \item `response`: Name of response variable (character)
#'     \item `true_params`: True parameter values (named numeric, may be NULL)
#'     \item Additional model-specific data
#'   }
#' @param fit_spec A list or single-row data.frame from fit_grid containing:
#'   \itemize{
#'     \item Model specification parameters
#'     \item Prior specifications
#'     \item Algorithm settings
#'   }
#' @param seed Integer seed for reproducibility
#' @param task_ctx A list with simulation context:
#'   \itemize{
#'     \item `task_id`: Unique task identifier
#'     \item `data_idx`: Index of data generation settings
#'     \item `fit_idx`: Index of fitting settings
#'     \item `rep_idx`: Replication index
#'   }
#'
#' @return A `bayesim_fit_result` S3 object containing:
#'   \itemize{
#'     \item success: TRUE/FALSE
#'     \item fit: backend object
#'     \item draws: matrix with colnames or NULL
#'     \item diagnostics: named list
#'     \item timing: list with total, warmup, sample
#'     \item warnings: character vector
#'     \item error: NULL or condition
#'   }
#' @export
fit <- S7::new_generic(
  "fit",
  "fitter",
  function(fitter, data_bundle, fit_spec, seed, task_ctx) {
    S7::S7_dispatch()
  }
)

#' @title Extract Posterior Draws
#' @description
#' Extract posterior draws from a fitted model.
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit()]
#' @param variables Character vector of variable names to extract.
#'   If NULL, extracts all variables.
#'
#' @return A matrix with dimensions S x P where:
#'   \itemize{
#'     \item S = number of posterior draws (iterations x chains)
#'     \item P = number of parameters/variables
#'     \item Column names match variable names
#'   }
#' @export
extract_draws <- S7::new_generic(
  "extract_draws",
  "fitter",
  function(fitter, fit_result, variables = NULL) {
    S7::S7_dispatch()
  }
)

#' @title Generate Predictions
#' @description
#' Generate predictions from a fitted model.
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit()]
#' @param newdata Data frame with new observations for prediction.
#'   If NULL, predictions are generated for the original training data.
#' @param seed Optional integer seed for reproducible predictions
#'
#' @return A list containing:
#'   \itemize{
#'     \item `predicted_mean`: Vector of mean predictions (N)
#'     \item `predicted_samples`: Matrix of posterior predictive samples (N x S)
#'     \item `predicted_sd`: Vector of prediction standard deviations (N)
#'     \item Additional fitter-specific outputs
#'   }
#' @export
predict_fit <- S7::new_generic(
  "predict_fit",
  "fitter",
  function(fitter, fit_result, newdata = NULL, seed = NULL) {
    S7::S7_dispatch()
  }
)

#' @title Compute Pointwise Log-Likelihood
#' @description
#' Compute pointwise log-likelihood values.
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit()]
#' @param newdata Data frame with observations. If NULL, uses training data.
#'
#' @return A matrix with dimensions N x S where:
#'   \itemize{
#'     \item N = number of observations
#'     \item S = number of posterior draws
#'     \item Entry (i, s) is log p(y_i | parameters_s)
#'   }
#' @export
log_lik <- S7::new_generic(
  "log_lik",
  "fitter",
  function(fitter, fit_result, newdata = NULL) {
    S7::S7_dispatch()
  }
)

#' @title Compute LOO-CV
#' @description
#' Compute leave-one-out cross-validation using Pareto-smoothed importance
#' sampling (PSIS-LOO).
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit()]
#'
#' @return A list containing:
#'   \itemize{
#'     \item `elpd`: Expected log predictive density (scalar)
#'     \item `p_loo`: Effective number of parameters (scalar)
#'     \item `elpd_se`: Standard error of ELPD (scalar)
#'     \item `pareto_k`: Pareto k diagnostic values (vector of length N)
#'     \item Additional loo-specific diagnostics
#'   }
#' @export
loo <- S7::new_generic("loo", "fitter", function(fitter, fit_result) {
  S7::S7_dispatch()
})

#' @title Extract Fit Diagnostics
#' @description
#' Extract convergence and fit diagnostics.
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit()]
#'
#' @return A named list of scalar diagnostic values, which may include:
#'   \itemize{
#'     \item `rhat_max`: Maximum R-hat statistic across parameters
#'     \item `ess_bulk`: Minimum bulk effective sample size
#'     \item `ess_tail`: Minimum tail effective sample size
#'     \item `divergent`: Number of divergent transitions
#'     \item `max_treedepth`: Number of max treedepth warnings
#'     \item Fitter-specific diagnostics
#'   }
#' @export
diagnostics <- S7::new_generic(
  "diagnostics",
  "fitter",
  function(fitter, fit_result) {
    S7::S7_dispatch()
  }
)


# =============================================================================
# Mock Fitter for Testing
# =============================================================================

#' @title Mock Fitter for Testing
#' @description
#' A simple implementation of the Fitter class for testing purposes.
#' This fitter returns predetermined values and does not perform actual
#' Bayesian inference.
#'
#' @return An S7 class object representing a MockFitter
#' @export
#' @seealso [Fitter] for the abstract base class
MockFitter <- S7::new_class(
  "MockFitter",
  parent = Fitter,
  properties = list(
    name = S7::new_property(S7::class_character, default = "mock"),
    supports_predictions = S7::new_property(S7::class_logical, default = TRUE),
    supports_log_lik = S7::new_property(S7::class_logical, default = TRUE),
    supports_loo = S7::new_property(S7::class_logical, default = TRUE),
    # Test-specific properties
    n_draws = S7::new_property(S7::class_integer, default = 100L),
    n_chains = S7::new_property(S7::class_integer, default = 4L)
  )
)

# =============================================================================
# MockFitter method implementations
# =============================================================================

S7::method(fit, MockFitter) <- function(
  fitter,
  data_bundle,
  fit_spec,
  seed,
  task_ctx
) {
  # Validate input
  if (is.null(data_bundle$train)) {
    rlang::abort(
      "data_bundle$train is required but got NULL",
      class = "bayesim_contract_error"
    )
  }

  n_obs <- nrow(data_bundle$train)
  timing_val <- runif(1, 0.5, 2.0) # Mock timing

  # Create mock fit object (stored in fit_result$fit)
  mock_fit <- list(
    fitter = "mock",
    data_bundle = data_bundle,
    fit_spec = fit_spec,
    task_ctx = task_ctx,
    n_obs = n_obs,
    seed = seed,
    timestamp = Sys.time()
  )

  new_fit_result(
    success = TRUE,
    fit = mock_fit,
    draws = extract_draws(
      fitter,
      list(seed = seed, data_bundle = data_bundle, n_obs = n_obs),
      variables = NULL
    ),
    diagnostics = list(
      rhat_max = 1.01,
      ess_bulk = 400,
      ess_tail = 350,
      divergent = 0L,
      max_treedepth = 0L
    ),
    timing = list(
      total = timing_val,
      warmup = timing_val / 2,
      sample = timing_val / 2
    ),
    warnings = character(),
    error = NULL
  )
}

S7::method(extract_draws, MockFitter) <- function(
  fitter,
  fit_result,
  variables = NULL
) {
  # Extract seed from fit object (fit_result$fit for proper structure,
  # or directly from fit_result for internal mock calls)
  fit_obj <- if (inherits(fit_result, "bayesim_fit_result")) {
    fit_result$fit
  } else {
    fit_result
  }
  seed <- fit_obj$seed %||% 12345L

  # Generate deterministic mock draws
  set.seed(seed)
  n_draws <- fitter@n_draws * fitter@n_chains

  # Base parameters
  param_names <- c("intercept", "beta", "sigma")
  if (!is.null(variables)) {
    param_names <- intersect(param_names, variables)
  }
  n_params <- length(param_names)

  # Create mock draws matrix
  draws <- matrix(
    rnorm(
      n_draws * n_params,
      mean = c(0, 1, 0.5)[seq_len(n_params)],
      sd = 0.1
    ),
    nrow = n_draws,
    ncol = n_params,
    dimnames = list(NULL, param_names)
  )

  # Ensure sigma is positive
  if ("sigma" %in% param_names) {
    draws[, "sigma"] <- abs(draws[, "sigma"])
  }

  draws
}

S7::method(predict_fit, MockFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL,
  seed = NULL
) {
  # Access data from fit_result$fit for proper structure
  fit_obj <- fit_result$fit
  data <- newdata %||% fit_obj$data_bundle$train
  n_obs <- nrow(data)

  # Use provided seed or fall back to fit seed
  rng_seed <- seed %||% fit_obj$seed %||% 12345L
  set.seed(rng_seed)

  # Get draws for prediction
  draws <- extract_draws(fitter, fit_result)
  n_draws <- nrow(draws)

  # Simple linear prediction: intercept + beta * 1 (mock feature)
  predicted_mean <- draws[, "intercept"] + draws[, "beta"]

  # Generate posterior predictive samples
  predicted_samples <- matrix(
    rnorm(n_obs * n_draws, mean = predicted_mean, sd = draws[, "sigma"]),
    nrow = n_obs,
    ncol = n_draws
  )

  list(
    predicted_mean = rowMeans(predicted_samples),
    predicted_samples = predicted_samples,
    predicted_sd = apply(predicted_samples, 1, sd)
  )
}

S7::method(log_lik, MockFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL
) {
  # Access data from fit_result$fit for proper structure
  fit_obj <- fit_result$fit
  data <- newdata %||% fit_obj$data_bundle$train
  n_obs <- nrow(data)

  set.seed(fit_obj$seed %||% 12345L)

  # Get draws
  draws <- extract_draws(fitter, fit_result)
  n_draws <- nrow(draws)

  # Simple normal log-likelihood
  # ll[i, s] = dnorm(y[i], mu[s], sigma[s], log = TRUE)
  ll <- matrix(
    rnorm(n_obs * n_draws, mean = -2, sd = 0.5), # Mock log-lik values
    nrow = n_obs,
    ncol = n_draws
  )

  ll
}

S7::method(loo, MockFitter) <- function(fitter, fit_result) {
  # Access data from fit_result$fit for proper structure
  fit_obj <- fit_result$fit
  n_obs <- fit_obj$n_obs
  set.seed(fit_obj$seed %||% 12345L)

  list(
    elpd = -n_obs * 2,
    p_loo = 3,
    elpd_se = sqrt(n_obs) * 0.5,
    pareto_k = runif(n_obs, -0.5, 0.5)
  )
}

S7::method(diagnostics, MockFitter) <- function(fitter, fit_result) {
  # Return diagnostics from fit_result if available, otherwise mock values
  if (!is.null(fit_result$diagnostics) && length(fit_result$diagnostics) > 0) {
    fit_result$diagnostics
  } else {
    list(
      rhat_max = 1.01,
      ess_bulk = 400,
      ess_tail = 350,
      divergent = 0L,
      max_treedepth = 0L
    )
  }
}


#' Null coalescing operator
#' @param x Value to check
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) if (is.null(x)) y else x
