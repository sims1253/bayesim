#' @title Abstract Fitter Class
#' @description
#' Abstract base class for Bayesian model fitters in bayesim.
#'
#' The Fitter class defines the interface that all model fitters must implement.
#' It provides a consistent API for fitting Bayesian models, extracting posterior
#' draws, generating predictions, computing log-likelihoods, and performing
#' model diagnostics.
#'
#' @param name Character string identifying the fitter (e.g., "stan", "brms")
#' @param supports_predictions Logical indicating if the fitter supports predictions
#' @param supports_log_lik Logical indicating if the fitter supports log-likelihood computation
#' @param supports_loo Logical indicating if the fitter supports LOO-CV
#'
#' @section Methods:
#' The following S7 generics must be implemented by subclasses:
#' \describe{
#'   \item{`fit_model(fitter, data_bundle, fit_spec, seed, task_ctx)`}{Main fitting method}
#'   \item{`extract_draws(fitter, fit_result, variables = NULL)`}{Extract posterior draws}
#'   \item{`predict_fit(fitter, fit_result, newdata = NULL, seed = NULL)`}{Generate predictions}
#'   \item{`log_lik_matrix(fitter, fit_result, newdata = NULL)`}{Pointwise log-likelihood}
#'   \item{`loo_fit(fitter, fit_result)`}{LOO-CV computation}
#'   \item{`fit_diagnostics(fitter, fit_result)`}{Extract fit diagnostics}
#' }
#'
#' @section Creating Custom Fitters:
#' To create a custom fitter, extend this class and implement methods for the
#' S7 generics: `fit_model()`, `extract_draws()`, `predict_fit()`,
#' `log_lik_matrix()`, `loo_fit()`, `fit_diagnostics()`. All matrices follow
#' the draws-by-observations (S x N) orientation; see
#' `vignette("custom-fitters")` for the full contract.
#'
#' @return An S7 class object representing the abstract Fitter
#' @export
#' @seealso [LinearRegressionFitter], [BrmsFitter], and [CmdStanFitter] for
#'   the built-in implementations, and [validate_fitter()] to check a custom
#'   fitter against the contract
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
fit_model <- S7::new_generic(
  "fit_model",
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
#' @param fit_result A `bayesim_fit_result` object from [fit_model()]
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
#' @param fit_result A `bayesim_fit_result` object from [fit_model()]
#' @param newdata Data frame with new observations for prediction.
#'   If NULL, predictions are generated for the original training data.
#' @param seed Optional integer seed for reproducible predictions
#'
#' @return A list containing:
#'   \itemize{
#'     \item `predicted_mean`: Vector of mean predictions (length N)
#'     \item `predicted_samples`: Matrix of posterior predictive samples
#'       (S x N; draws as rows, observations as columns)
#'     \item `predicted_sd`: Vector of prediction standard deviations (length N)
#'     \item Additional fitter-specific outputs
#'   }
#'
#'   `predicted_samples` follows the same orientation convention as
#'   [log_lik_matrix()] and [predict_epred()]: all matrices are draws x observations
#'   (S rows, N columns).
#' @export
predict_fit <- S7::new_generic(
  "predict_fit",
  "fitter",
  function(fitter, fit_result, newdata = NULL, seed = NULL) {
    S7::S7_dispatch()
  }
)

#' @title Compute Posterior Expectation Predictions
#' @description
#' Compute posterior draws of the expected value of the response distribution
#' (mu, without observation noise). This is the `epred` quantity used by brms'
#' `loo_R2()` and bayesim's `r2_loo_metric()` (F3). It must NOT include
#' observation-level noise — only the model's conditional mean.
#'
#' Fitters that cannot provide expectation predictions should return `NULL`; the
#' `r2_loo` metric then degrades to NA. The default Fitter method returns NULL.
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit_model()]
#' @param newdata Data frame with observations. If NULL, uses training data.
#'
#' @return A matrix with dimensions S x N (draws x observations), or NULL if
#'   not supported by the fitter.
#' @export
predict_epred <- S7::new_generic(
  "predict_epred",
  "fitter",
  function(fitter, fit_result, newdata = NULL) {
    S7::S7_dispatch()
  }
)

# Default: fitters that don't override return NULL (metric degrades to NA).
S7::method(predict_epred, S7::class_any) <- function(
  fitter,
  fit_result,
  newdata = NULL
) {
  NULL
}

#' @title Compute Pointwise Log-Likelihood
#' @description
#' Compute pointwise log-likelihood values. Named `log_lik_matrix` (rather than
#' `log_lik`) so that bayesim does not mask [brms::log_lik] or the rstantools
#' `log_lik` generic for users who load bayesim alongside brms.
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit_model()]
#' @param newdata Data frame with observations. If NULL, uses training data.
#'
#' @return A matrix with dimensions S x N where:
#'   \itemize{
#'     \item S = number of posterior draws (rows)
#'     \item N = number of observations (columns)
#'     \item Entry (s, i) is log p(y_i | parameters_s)
#'   }
#'   This is the brms/loo convention (draws as rows), matching
#'   `brms::log_lik` and what `loo::psis`/`loo::E_loo`/`loo::relative_eff`
#'   expect.
#' @export
log_lik_matrix <- S7::new_generic(
  "log_lik_matrix",
  "fitter",
  function(fitter, fit_result, newdata = NULL) {
    S7::S7_dispatch()
  }
)

#' @title Compute LOO-CV
#' @description
#' Compute leave-one-out cross-validation using Pareto-smoothed importance
#' sampling (PSIS-LOO). Named `loo_fit` to avoid clashing with [loo::loo()].
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit_model()]
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
loo_fit <- S7::new_generic("loo_fit", "fitter", function(fitter, fit_result) {
  S7::S7_dispatch()
})

#' @title Extract Fit Diagnostics
#' @description
#' Extract convergence and fit diagnostics. Named `fit_diagnostics` (rather
#' than `diagnostics`) to avoid exporting a generic-noun that collides with
#' other packages.
#'
#' @param fitter An S7 Fitter object
#' @param fit_result A `bayesim_fit_result` object from [fit_model()]
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
fit_diagnostics <- S7::new_generic(
  "fit_diagnostics",
  "fitter",
  function(fitter, fit_result) {
    S7::S7_dispatch()
  }
)


# =============================================================================
# Mock Fitter for Testing
# =============================================================================

#' @title Mock Fitter for Testing (Internal)
#' @description
#' A simple implementation of the Fitter class for **testing purposes only**.
#' This fitter returns predetermined values and does not perform actual
#' Bayesian inference.
#'
#' **WARNING**: MockFitter ignores the `fit_spec` model specification.
#' Do NOT use for scientific simulation studies comparing different models.
#' For real inference, use [BrmsFitter] or create a custom fitter.
#'
#' For a lightweight no-Stan alternative, see the `LinearFitter` example
#' in `vignette("custom-fitters")`.
#'
#' @param name Character string identifying the fitter
#' @param supports_predictions Logical; whether predictions are supported
#' @param supports_log_lik Logical; whether log-likelihood is supported
#' @param supports_loo Logical; whether LOO-CV is supported
#' @param n_draws Integer; number of posterior draws to simulate
#' @param n_chains Integer; number of chains to simulate
#'
#' @return An S7 class object representing a MockFitter
#' @keywords internal
#' @seealso [Fitter] for the abstract base class, [BrmsFitter] and
#'   [LinearRegressionFitter] for real inference
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

S7::method(fit_model, MockFitter) <- function(
  fitter,
  data_bundle,
  fit_spec,
  seed,
  task_ctx
) {
  # Warn users this is for testing only
  if (!is.null(task_ctx$fit_idx) && task_ctx$fit_idx > 1) {
    cli::cli_warn(c(
      "MockFitter ignores fit_spec and is for testing only.",
      "i" = "Do not use for model comparison studies.",
      "i" = "Use BrmsFitter or a custom fitter for real inference."
    ))
  }

  # Validate input
  if (is.null(data_bundle$train)) {
    rlang::abort(
      "data_bundle$train is required but got NULL",
      class = "bayesim_contract_error"
    )
  }

  n_obs <- nrow(data_bundle$train)

  # Use with_seed to avoid polluting global RNG state
  timing_val <- withr::with_seed(seed %||% 12345L, {
    runif(1, 0.5, 2.0) # Mock timing
  })

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
    error = NULL,
    data_bundle = data_bundle
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
  data_bundle <- fit_obj$data_bundle

  # Generate deterministic mock draws without polluting global RNG
  withr::with_seed(seed, {
    n_draws <- fitter@n_draws * fitter@n_chains

    # Get parameter names with fallback hierarchy:
    # 1. vars_of_interest from data_bundle
    # 2. names of true_params from data_bundle
    # 3. hardcoded defaults
    default_params <- c("intercept", "beta", "sigma")

    param_names <- data_bundle$vars_of_interest
    if (is.null(param_names)) {
      param_names <- names(data_bundle$true_params)
    }
    if (is.null(param_names)) {
      param_names <- default_params
    }

    if (!is.null(variables)) {
      param_names <- intersect(param_names, variables)
    }
    n_params <- length(param_names)

    # Get true parameter values for centering draws
    true_params <- data_bundle$true_params %||% list()

    # Determine center values for each parameter
    # Use true_params if available, otherwise use sensible defaults:
    # - 0 for location parameters (intercept, beta, etc.)
    # - 1 for scale parameters (sigma, sd, scale, etc.)
    center_values <- vapply(
      param_names,
      function(param) {
        if (!is.null(true_params[[param]])) {
          true_params[[param]]
        } else if (grepl("sigma|sd|scale|var", param, ignore.case = TRUE)) {
          1 # Default for scale parameters
        } else {
          0 # Default for location parameters
        }
      },
      numeric(1)
    )

    # Create mock draws matrix centered around true/default values
    draws <- matrix(
      rnorm(n_draws * n_params, mean = center_values, sd = 0.1),
      nrow = n_draws,
      ncol = n_params,
      dimnames = list(NULL, param_names)
    )

    # Ensure scale parameters are positive
    scale_params <- grep(
      "sigma|sd|scale|var",
      param_names,
      ignore.case = TRUE,
      value = TRUE
    )
    for (param in scale_params) {
      draws[, param] <- abs(draws[, param])
    }

    draws
  })
}

S7::method(predict_fit, MockFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL,
  seed = NULL
) {
  # Access data from fit_result$data_bundle (preferred) or fit_obj$data_bundle (fallback)
  fit_obj <- fit_result$fit
  data_bundle <- fit_result$data_bundle %||% fit_obj$data_bundle
  data <- newdata %||% data_bundle$train
  n_obs <- nrow(data)

  # Use provided seed or fall back to fit seed
  rng_seed <- seed %||% fit_obj$seed %||% 12345L

  withr::with_seed(rng_seed, {
    # Get draws for prediction
    draws <- extract_draws(fitter, fit_result)
    n_draws <- nrow(draws)
    param_names <- colnames(draws)

    # Find parameter names dynamically using pattern matching
    find_param <- function(pattern, params) {
      matches <- grep(pattern, params, ignore.case = TRUE, value = TRUE)
      if (length(matches) > 0) matches[1] else NULL
    }

    intercept_param <- find_param("intercept|alpha|a_", param_names)
    slope_param <- find_param("slope|beta|b_", param_names)
    scale_param <- find_param("sigma|sd|scale", param_names)

    # Determine predictor values (x)
    # Get predictor columns (exclude response variable if known)
    response_var <- data_bundle$response
    predictor_cols <- setdiff(names(data), response_var)

    if (length(predictor_cols) > 0) {
      # Use the first predictor column
      x <- data[[predictor_cols[1]]]
    } else {
      # No predictor available, use constant 1
      x <- 1
    }

    # Calculate per-draw, per-observation mean as an S x N matrix
    # (draws as rows, observations as columns). draws is S x P and x is
    # length N, so draws[, p] * x recycles column-wise into S x N.
    if (!is.null(intercept_param) && !is.null(slope_param)) {
      # Linear prediction: intercept + slope * x  -> S x N
      intercept_vec <- draws[, intercept_param] # length S
      slope_vec <- draws[, slope_param] # length S
      pred_mean_matrix <- outer(intercept_vec, rep(1, n_obs)) +
        outer(slope_vec, x)
    } else if (!is.null(intercept_param)) {
      # Only intercept available: replicate intercept across observations.
      intercept_vec <- draws[, intercept_param]
      pred_mean_matrix <- matrix(
        rep(intercept_vec, n_obs),
        nrow = n_draws,
        ncol = n_obs
      )
    } else if (!is.null(slope_param)) {
      # Only slope available (no intercept): slope * x  -> S x N
      slope_vec <- draws[, slope_param]
      pred_mean_matrix <- outer(slope_vec, x)
    } else {
      # No location parameters found, use 0 as fallback (S x N of zeros).
      pred_mean_matrix <- matrix(0, nrow = n_draws, ncol = n_obs)
    }

    # Get scale values for predictions (length-S vector, one per draw)
    if (!is.null(scale_param)) {
      scale_values <- draws[, scale_param]
    } else {
      # Default scale if not found
      scale_values <- rep(1, n_draws)
    }

    # Build the S x N scale matrix by recycling the length-S scale_values
    # across the N columns (each draw has one scale applied to all obs).
    scale_matrix <- matrix(
      rep(scale_values, n_obs),
      nrow = n_draws,
      ncol = n_obs
    )

    # Generate posterior predictive samples as S x N
    # (draws as rows, observations as columns).
    predicted_samples <- matrix(
      rnorm(n_draws * n_obs, mean = pred_mean_matrix, sd = scale_matrix),
      nrow = n_draws,
      ncol = n_obs
    )

    list(
      predicted_mean = colMeans(predicted_samples),
      predicted_samples = predicted_samples,
      predicted_sd = apply(predicted_samples, 2, sd)
    )
  })
}

S7::method(log_lik_matrix, MockFitter) <- function(
  fitter,
  fit_result,
  newdata = NULL
) {
  # Access data from fit_result$data_bundle (preferred) or fit_obj$data_bundle (fallback)
  fit_obj <- fit_result$fit
  data_bundle <- fit_result$data_bundle %||% fit_obj$data_bundle
  data <- newdata %||% data_bundle$train
  n_obs <- nrow(data)

  withr::with_seed(fit_obj$seed %||% 12345L, {
    # Get draws
    draws <- extract_draws(fitter, fit_result)
    n_draws <- nrow(draws)

    # Simple normal log-likelihood.
    # Convention: S x N matrix (draws x observations), matching brms/loo.
    ll <- matrix(
      rnorm(n_obs * n_draws, mean = -2, sd = 0.5), # Mock log-lik values
      nrow = n_draws,
      ncol = n_obs
    )

    ll
  })
}

S7::method(loo_fit, MockFitter) <- function(fitter, fit_result) {
  if (is.null(fit_result) || is.null(fit_result$fit)) {
    return(list(
      elpd = NA_real_,
      p_loo = NA_real_,
      elpd_se = NA_real_,
      pareto_k = numeric()
    ))
  }

  # Access data from fit_result$fit for proper structure
  fit_obj <- fit_result$fit
  n_obs <- fit_obj$n_obs

  withr::with_seed(fit_obj$seed %||% 12345L, {
    list(
      elpd = -n_obs * 2,
      p_loo = 3,
      elpd_se = sqrt(n_obs) * 0.5,
      pareto_k = runif(n_obs, -0.5, 0.5)
    )
  })
}

S7::method(fit_diagnostics, MockFitter) <- function(fitter, fit_result) {
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


#' @title Validate a Fitter Object
#' @description
#' Validates that a fitter object correctly implements the bayesim Fitter interface.
#' This function checks that the object is a valid S7 Fitter class instance with
#' all required properties and methods implemented.
#'
#' @param fitter An S7 Fitter object to validate
#' @param smoke_test Logical, if TRUE run a quick fit test with sample data to
#'   verify that methods work correctly end-to-end
#' @param verbose Logical, if TRUE print progress messages during validation
#'
#' @return The validated fitter object (invisibly) if valid, otherwise raises an error with details about what failed
#'
#' @details
#' The validation performs the following checks:
#'
#' **Property Checks:**
#' - Object is an S7 Fitter class
#' - `name` property exists and is character
#' - `supports_predictions` property exists and is logical
#' - `supports_log_lik` property exists and is logical
#' - `supports_loo` property exists and is logical
#'
#' **Method Checks:**
#' - `fit_model()` method is implemented
#' - `extract_draws()` method is implemented
#' - `predict_fit()` method is implemented
#' - `log_lik_matrix()` method is implemented
#' - `loo_fit()` method is implemented
#' - `fit_diagnostics()` method is implemented
#'
#' **Smoke Test (when smoke_test = TRUE):**
#' - Creates simple lm-like test data
#' - Calls `fit_model()` and verifies `bayesim_fit_result` structure
#' - Calls `extract_draws()` and verifies matrix with colnames
#' - If `supports_predictions`, calls `predict_fit()` and verifies output
#' - If `supports_log_lik`, calls `log_lik_matrix()` and verifies matrix output
#' - Calls `fit_diagnostics()` and verifies list output
#'
#' @export
#' @seealso [Fitter], [MockFitter], [validate_metric()]
#'
#' @examples
#' # Validate a built-in fitter (basic check only)
#' validate_fitter(LinearRegressionFitter())
#'
#' # Full validation with an end-to-end smoke test
#' validate_fitter(LinearRegressionFitter(), smoke_test = TRUE)
#'
#' \dontrun{
#' # Use in your own fitter's tests
#' my_fitter <- MyCustomFitter()
#' validate_fitter(my_fitter, smoke_test = TRUE, verbose = TRUE)
#' }
validate_fitter <- function(fitter, smoke_test = FALSE, verbose = FALSE) {
  # Helper function for conditional messages
  msg <- function(...) {
    if (verbose) {
      message(...)
    }
  }

  # ===========================================================================
  # Check 1: S7 Fitter class
  # ===========================================================================
  msg("Checking S7 Fitter class...")
  if (!S7::S7_inherits(fitter, Fitter)) {
    rlang::abort(
      c(
        "Object is not an S7 Fitter class",
        i = "Object class: ",
        paste(class(fitter), collapse = ", "),
        i = "Use S7::new_class() with parent = Fitter to create a fitter"
      ),
      class = "bayesim_validation_error"
    )
  }
  msg("  [OK] Object is an S7 Fitter class")

  # ===========================================================================
  # Check 2: Required properties
  # ===========================================================================
  msg("Checking required properties...")

  required_props <- c(
    "name",
    "supports_predictions",
    "supports_log_lik",
    "supports_loo"
  )

  for (prop in required_props) {
    # Use S7 property accessor via @ operator
    prop_value <- tryCatch(
      {
        eval(substitute(fitter@p, list(p = prop)))
      },
      error = function(e) NULL
    )

    if (is.null(prop_value)) {
      rlang::abort(
        c(
          paste0("Missing required property: '", prop, "'"),
          i = "All Fitter subclasses must have: name, supports_predictions, supports_log_lik, supports_loo"
        ),
        class = "bayesim_validation_error"
      )
    }
    msg("  [OK] Property '", prop, "' exists: ", toString(prop_value))
  }

  # ===========================================================================
  # Check 3: Required S7 methods
  # ===========================================================================
  msg("Checking required S7 methods...")

  required_methods <- c(
    "fit_model",
    "extract_draws",
    "predict_fit",
    "log_lik_matrix",
    "loo_fit",
    "fit_diagnostics"
  )
  fitter_class <- S7::S7_class(fitter)

  for (method_name in required_methods) {
    # Check if the method is defined for this fitter's class
    method_impl <- tryCatch(
      S7::method(get(method_name), fitter_class),
      error = function(e) NULL
    )

    if (is.null(method_impl)) {
      rlang::abort(
        c(
          paste0("Missing required method: '", method_name, "()'"),
          i = paste0(
            "Use S7::method(",
            method_name,
            ", YourFitterClass) <- function(fitter, ...) { ... }"
          )
        ),
        class = "bayesim_validation_error"
      )
    }
    msg("  [OK] Method '", method_name, "()' is implemented")
  }

  # ===========================================================================
  # Check 4: Smoke test (optional)
  # ===========================================================================
  if (smoke_test) {
    msg("Running smoke test with sample data...")

    # Create simple test data (lm-like)
    set.seed(12345)
    n <- 20

    test_data_bundle <- list(
      train = data.frame(
        y = rnorm(n),
        x = rnorm(n)
      ),
      test = NULL,
      response = "y",
      true_params = c(intercept = 0, beta = 0, sigma = 1),
      vars_of_interest = c("intercept", "beta", "sigma"),
      meta = list()
    )

    test_fit_spec <- data.frame(model = "test")
    test_task_ctx <- list(
      task_id = "smoke_test",
      data_idx = 1L,
      fit_idx = 1L,
      rep_idx = 1L
    )

    # Test fit_model()
    msg("  Testing fit_model()...")
    fit_result <- tryCatch(
      fit_model(
        fitter,
        test_data_bundle,
        test_fit_spec,
        seed = 12345L,
        test_task_ctx
      ),
      error = function(e) {
        rlang::abort(
          c(
            "fit_model() method failed during smoke test",
            x = conditionMessage(e)
          ),
          class = "bayesim_validation_error"
        )
      }
    )

    # Validate fit_result structure
    if (!inherits(fit_result, "bayesim_fit_result")) {
      rlang::abort(
        c(
          "fit_model() did not return a bayesim_fit_result object",
          i = paste0(
            "Returned class: ",
            paste(class(fit_result), collapse = ", ")
          )
        ),
        class = "bayesim_validation_error"
      )
    }

    if (!isTRUE(fit_result$success)) {
      rlang::abort(
        c(
          "fit_model() returned unsuccessful result during smoke test",
          i = if (!is.null(fit_result$error)) {
            conditionMessage(fit_result$error)
          } else {
            "No error message"
          }
        ),
        class = "bayesim_validation_error"
      )
    }
    msg("    [OK] fit_model() returns valid bayesim_fit_result")

    # Test extract_draws()
    msg("  Testing extract_draws()...")
    draws <- tryCatch(
      extract_draws(fitter, fit_result),
      error = function(e) {
        rlang::abort(
          c(
            "extract_draws() method failed during smoke test",
            x = conditionMessage(e)
          ),
          class = "bayesim_validation_error"
        )
      }
    )

    if (!is.matrix(draws)) {
      rlang::abort(
        c(
          "extract_draws() did not return a matrix",
          i = paste0("Returned class: ", paste(class(draws), collapse = ", "))
        ),
        class = "bayesim_validation_error"
      )
    }

    if (is.null(colnames(draws))) {
      rlang::abort(
        c(
          "extract_draws() returned matrix without column names",
          i = "Matrix columns should be named after parameters"
        ),
        class = "bayesim_validation_error"
      )
    }
    msg("    [OK] extract_draws() returns matrix with colnames")

    # Test predict_fit() if supported
    if (isTRUE(fitter@supports_predictions)) {
      msg("  Testing predict_fit()...")
      preds <- tryCatch(
        predict_fit(fitter, fit_result),
        error = function(e) {
          rlang::abort(
            c(
              "predict_fit() method failed during smoke test",
              x = conditionMessage(e)
            ),
            class = "bayesim_validation_error"
          )
        }
      )

      if (is.null(preds)) {
        rlang::abort(
          c(
            "predict_fit() returned NULL but supports_predictions is TRUE",
            i = "Return a list with predicted_mean, predicted_samples, and predicted_sd"
          ),
          class = "bayesim_validation_error"
        )
      }

      if (
        is.null(preds$predicted_mean) ||
          is.null(preds$predicted_samples) ||
          is.null(preds$predicted_sd)
      ) {
        rlang::abort(
          c(
            "predict_fit() result missing required elements",
            i = "Must contain: predicted_mean, predicted_samples, predicted_sd"
          ),
          class = "bayesim_validation_error"
        )
      }

      if (length(preds$predicted_mean) != n) {
        rlang::abort(
          c(
            "predict_fit() returned wrong number of predictions",
            i = paste0("Expected ", n, ", got ", length(preds$predicted_mean))
          ),
          class = "bayesim_validation_error"
        )
      }

      # Orientation check: predicted_samples must be S x N (draws x
      # observations). We verify observations-as-columns (N columns). The row
      # count S is fitter-specific, so only the column count is enforced; this
      # is robust to fitters whose n_draws happens to equal n_obs by coincidence.
      if (is.matrix(preds$predicted_samples)) {
        if (ncol(preds$predicted_samples) != n) {
          rlang::abort(
            c(
              "predict_fit() returned predicted_samples with the wrong orientation",
              i = paste0(
                "predicted_samples must be S x N (draws x observations); ",
                "expected N columns, got ",
                ncol(preds$predicted_samples)
              )
            ),
            class = "bayesim_validation_error"
          )
        }
      }
      msg("    [OK] predict_fit() returns valid predictions")
    } else {
      msg("  [SKIP] predict_fit() (supports_predictions is FALSE)")
    }

    # Test log_lik_matrix() if supported.
    # Convention: log_lik_matrix() must return an S x N matrix (draws x observations),
    # matching the brms/loo orientation. We check that the number of columns
    # equals n_obs (one column per observation). The S rows are posterior
    # draws and their count is fitter-specific, so the row count is not
    # checked here; this is robust to fitters whose n_draws happens to equal
    # n_obs only by coincidence.
    if (isTRUE(fitter@supports_log_lik)) {
      msg("  Testing log_lik_matrix()...")
      ll <- tryCatch(
        log_lik_matrix(fitter, fit_result),
        error = function(e) {
          rlang::abort(
            c(
              "log_lik_matrix() method failed during smoke test",
              x = conditionMessage(e)
            ),
            class = "bayesim_validation_error"
          )
        }
      )

      if (is.null(ll)) {
        rlang::abort(
          c(
            "log_lik_matrix() returned NULL but supports_log_lik is TRUE",
            i = "Return an S x N matrix of pointwise log-likelihoods (draws x observations)"
          ),
          class = "bayesim_validation_error"
        )
      }

      if (!is.matrix(ll)) {
        rlang::abort(
          c(
            "log_lik_matrix() did not return a matrix",
            i = paste0("Returned class: ", paste(class(ll), collapse = ", "))
          ),
          class = "bayesim_validation_error"
        )
      }

      if (ncol(ll) != n) {
        rlang::abort(
          c(
            "log_lik_matrix() returned wrong number of columns",
            i = paste0(
              "Expected N columns (one per observation), got ",
              ncol(ll),
              " (expected ",
              n,
              ")"
            )
          ),
          class = "bayesim_validation_error"
        )
      }
      msg("    [OK] log_lik_matrix() returns valid matrix")
    } else {
      msg("  [SKIP] log_lik_matrix() (supports_log_lik is FALSE)")
    }

    # Test fit_diagnostics()
    msg("  Testing fit_diagnostics()...")
    diag <- tryCatch(
      fit_diagnostics(fitter, fit_result),
      error = function(e) {
        rlang::abort(
          c(
            "fit_diagnostics() method failed during smoke test",
            x = conditionMessage(e)
          ),
          class = "bayesim_validation_error"
        )
      }
    )

    if (!is.list(diag)) {
      rlang::abort(
        c(
          "fit_diagnostics() did not return a list",
          i = paste0("Returned class: ", paste(class(diag), collapse = ", "))
        ),
        class = "bayesim_validation_error"
      )
    }
    msg("    [OK] fit_diagnostics() returns list")

    msg("Smoke test completed successfully!")
  }

  msg("Validation passed!")
  invisible(fitter)
}
