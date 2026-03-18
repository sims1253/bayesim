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
#'   \item{`fit(fitter, data_bundle, fit_spec, seed, task_ctx)`}{Main fitting method}
#'   \item{`extract_draws(fitter, fit_result, variables = NULL)`}{Extract posterior draws}
#'   \item{`predict_fit(fitter, fit_result, newdata = NULL, seed = NULL)`}{Generate predictions}
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
#' @return An S7 class object representing a MockFitter
#' @export
#' @keywords internal
#' @seealso [Fitter] for the abstract base class, [BrmsFitter] for real inference
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
    stop(bayesim_contract_error("data_bundle$train is required but got NULL"))
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

    # Calculate predicted mean based on available parameters
    if (!is.null(intercept_param) && !is.null(slope_param)) {
      # Linear prediction: intercept + slope * x
      # Result is a vector of length n_draws (for scalar x) or matrix (for vector x)
      predicted_mean <- draws[, intercept_param] + draws[, slope_param] * x
    } else if (!is.null(intercept_param)) {
      # Only intercept available
      predicted_mean <- draws[, intercept_param]
    } else if (!is.null(slope_param)) {
      # Only slope available (no intercept)
      predicted_mean <- draws[, slope_param] * x
    } else {
      # No location parameters found, use 0 as fallback
      predicted_mean <- rep(0, n_draws)
    }

    # Get scale values for predictions
    if (!is.null(scale_param)) {
      scale_values <- draws[, scale_param]
    } else {
      # Default scale if not found
      scale_values <- rep(1, n_draws)
    }

    # Generate posterior predictive samples
    # predicted_mean should be n_draws-length for scalar x,
    # or n_obs x n_draws matrix for vector x
    if (is.matrix(predicted_mean)) {
      # Already a matrix (n_obs x n_draws)
      pred_mean_matrix <- predicted_mean
    } else {
      # Replicate to n_obs x n_draws matrix
      pred_mean_matrix <- matrix(
        rep(predicted_mean, each = n_obs),
        nrow = n_obs,
        ncol = n_draws
      )
    }

    predicted_samples <- matrix(
      rnorm(n_obs * n_draws, mean = pred_mean_matrix, sd = scale_values),
      nrow = n_obs,
      ncol = n_draws
    )

    list(
      predicted_mean = rowMeans(predicted_samples),
      predicted_samples = predicted_samples,
      predicted_sd = apply(predicted_samples, 1, sd)
    )
  })
}

S7::method(log_lik, MockFitter) <- function(
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

    # Simple normal log-likelihood
    ll <- matrix(
      rnorm(n_obs * n_draws, mean = -2, sd = 0.5), # Mock log-lik values
      nrow = n_obs,
      ncol = n_draws
    )

    ll
  })
}

S7::method(loo, MockFitter) <- function(fitter, fit_result) {
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
#' - `fit()` method is implemented
#' - `extract_draws()` method is implemented
#' - `predict_fit()` method is implemented
#' - `log_lik()` method is implemented
#' - `loo()` method is implemented
#' - `diagnostics()` method is implemented
#'
#' **Smoke Test (when smoke_test = TRUE):**
#' - Creates simple lm-like test data
#' - Calls `fit()` and verifies `bayesim_fit_result` structure
#' - Calls `extract_draws()` and verifies matrix with colnames
#' - If `supports_predictions`, calls `predict_fit()` and verifies output
#' - If `supports_log_lik`, calls `log_lik()` and verifies matrix output
#' - Calls `diagnostics()` and verifies list output
#'
#' @export
#' @seealso [Fitter], [MockFitter]
#'
#' @examples
#' # Validate the built-in MockFitter (basic check only)
#' validate_fitter(MockFitter())
#'
#' # Full validation with smoke test
#' validate_fitter(MockFitter(), smoke_test = TRUE, verbose = TRUE)
#'
#' # Use in your fitter tests
#' \dontrun{
#' my_fitter <- MyCustomFitter()
#' validate_fitter(my_fitter, smoke_test = TRUE)
#' cat("Fitter is valid!\n")
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
    stop(bayesim_validation_error(
      paste0(
        "Object is not an S7 Fitter class. ",
        "Object class: ", paste(class(fitter), collapse = ", "), ". ",
        "Use S7::new_class() with parent = Fitter to create a fitter."
      )
    ))
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
      stop(bayesim_validation_error(
        paste0("Missing required property: '", prop, "'. ",
              "All Fitter subclasses must have: name, supports_predictions, supports_log_lik, supports_loo")
      ))
    }
    msg("  [OK] Property '", prop, "' exists: ", toString(prop_value))
  }

  # ===========================================================================
  # Check 3: Required S7 methods
  # ===========================================================================
  msg("Checking required S7 methods...")

  required_methods <- c(
    "fit",
    "extract_draws",
    "predict_fit",
    "log_lik",
    "loo",
    "diagnostics"
  )
  fitter_class <- S7::S7_class(fitter)

  for (method_name in required_methods) {
    # Check if the method is defined for this fitter's class
    method_impl <- tryCatch(
      S7::method(get(method_name), fitter_class),
      error = function(e) NULL
    )

    if (is.null(method_impl)) {
      stop(bayesim_validation_error(
        paste0("Missing required method: '", method_name, "()'. ",
              "Use S7::method(", method_name, ", YourFitterClass) <- function(fitter, ...) { ... }")
      ))
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
      references = NULL,
      meta = list()
    )

    test_fit_spec <- data.frame(model = "test")
    test_task_ctx <- list(
      task_id = "smoke_test",
      data_idx = 1L,
      fit_idx = 1L,
      rep_idx = 1L
    )

    # Test fit()
    msg("  Testing fit()...")
    fit_result <- tryCatch(
      fit(
        fitter,
        test_data_bundle,
        test_fit_spec,
        seed = 12345L,
        test_task_ctx
      ),
      error = function(e) {
        stop(bayesim_validation_error(
          paste0("fit() method failed during smoke test: ", conditionMessage(e))
        ))
      }
    )

    # Validate fit_result structure
    if (!inherits(fit_result, "bayesim_fit_result")) {
      stop(bayesim_validation_error(
        paste0("fit() did not return a bayesim_fit_result object. ",
              "Returned class: ", paste(class(fit_result), collapse = ", "))
      ))
    }

    if (!isTRUE(fit_result$success)) {
      stop(bayesim_validation_error(
        paste0("fit() returned unsuccessful result during smoke test: ",
              if (!is.null(fit_result$error)) conditionMessage(fit_result$error) else "No error message")
      ))
    }
    msg("    [OK] fit() returns valid bayesim_fit_result")

    # Test extract_draws()
    msg("  Testing extract_draws()...")
    draws <- tryCatch(
      extract_draws(fitter, fit_result),
      error = function(e) {
        stop(bayesim_validation_error(
          paste0("extract_draws() method failed during smoke test: ", conditionMessage(e))
        ))
      }
    )

    if (!is.matrix(draws)) {
      stop(bayesim_validation_error(
        paste0("extract_draws() did not return a matrix. Returned class: ", paste(class(draws), collapse = ", "))
      ))
    }

    if (is.null(colnames(draws))) {
      stop(bayesim_validation_error(
        "extract_draws() returned matrix without column names. Matrix columns should be named after parameters"
      ))
    }
    msg("    [OK] extract_draws() returns matrix with colnames")

    # Test predict_fit() if supported
    if (isTRUE(fitter@supports_predictions)) {
      msg("  Testing predict_fit()...")
      preds <- tryCatch(
        predict_fit(fitter, fit_result),
        error = function(e) {
          stop(bayesim_validation_error(
            paste0("predict_fit() method failed during smoke test: ", conditionMessage(e))
          ))
        }
      )

      if (is.null(preds)) {
        stop(bayesim_validation_error(
          "predict_fit() returned NULL but supports_predictions is TRUE. Return a list with predicted_mean, predicted_samples, and predicted_sd"
        ))
      }

      if (
        is.null(preds$predicted_mean) ||
          is.null(preds$predicted_samples) ||
          is.null(preds$predicted_sd)
      ) {
        stop(bayesim_validation_error(
          "predict_fit() result missing required elements. Must contain: predicted_mean, predicted_samples, predicted_sd"
        ))
      }

      if (length(preds$predicted_mean) != n) {
        stop(bayesim_validation_error(
          paste0("predict_fit() returned wrong number of predictions. Expected ", n, ", got ", length(preds$predicted_mean))
        ))
      }
      msg("    [OK] predict_fit() returns valid predictions")
    } else {
      msg("  [SKIP] predict_fit() (supports_predictions is FALSE)")
    }

    # Test log_lik() if supported
    if (isTRUE(fitter@supports_log_lik)) {
      msg("  Testing log_lik()...")
      ll <- tryCatch(
        log_lik(fitter, fit_result),
        error = function(e) {
          stop(bayesim_validation_error(
            paste0("log_lik() method failed during smoke test: ", conditionMessage(e))
          ))
        }
      )

      if (is.null(ll)) {
        stop(bayesim_validation_error(
          "log_lik() returned NULL but supports_log_lik is TRUE. Return an N x S matrix of pointwise log-likelihoods"
        ))
      }

      if (!is.matrix(ll)) {
        stop(bayesim_validation_error(
          paste0("log_lik() did not return a matrix. Returned class: ", paste(class(ll), collapse = ", "))
        ))
      }

      if (nrow(ll) != n) {
        stop(bayesim_validation_error(
          paste0("log_lik() returned wrong number of rows. Expected ", n, " rows (one per observation), got ", nrow(ll))
        ))
      }
      msg("    [OK] log_lik() returns valid matrix")
    } else {
      msg("  [SKIP] log_lik() (supports_log_lik is FALSE)")
    }

    # Test diagnostics()
    msg("  Testing diagnostics()...")
    diag <- tryCatch(
      diagnostics(fitter, fit_result),
      error = function(e) {
        stop(bayesim_validation_error(
          paste0("diagnostics() method failed during smoke test: ", conditionMessage(e))
        ))
      }
    )

    if (!is.list(diag)) {
      stop(bayesim_validation_error(
        paste0("diagnostics() did not return a list. Returned class: ", paste(class(diag), collapse = ", "))
      ))
    }
    msg("    [OK] diagnostics() returns list")

    msg("Smoke test completed successfully!")
  }

  msg("Validation passed!")
  invisible(fitter)
}
