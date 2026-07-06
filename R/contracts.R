#' @title Contract Validators
#' @description Validators for data bundles, fit results, fitter interfaces,
#'   metric interfaces, and simulation configurations. These validators enforce
#'   the contracts between components in the bayesim package.
#' @name contracts
#' @keywords internal
NULL

# =============================================================================
# Data Bundle Validation
# =============================================================================

#' Validate a Data Bundle
#'
#' Validates that a data_bundle conforms to the required structure for use
#' in simulation tasks. The data_bundle is the output of a data_generator
#' function and contains all data-related objects needed for model fitting
#' and metric computation.
#'
#' @param data_bundle A list containing data and metadata for a simulation task.
#'
#' @return The input `data_bundle`, invisibly, if validation passes.
#'
#' @details
#' The data_bundle must have the following structure:
#' \itemize{
#'   \item `train`: A data.frame with at least 1 row (required)
#'   \item `test`: NULL or a data.frame (optional)
#'   \item `response`: A scalar character naming the response column in train
#'     (and test if not NULL)
#'   \item `true_params`: A named numeric vector where names exactly match
#'     `vars_of_interest`
#'   \item `vars_of_interest`: A non-empty character vector of unique names
#'   \item `meta`: Optional named list with scalar values only
#' }
#'
#' Validation rules:
#' \itemize{
#'   \item `train` must be a data.frame with nrow >= 1
#'   \item `test` must be NULL or a data.frame
#'   \item `response` must be a scalar character present in train (and test if not NULL)
#'   \item `true_params` must be a named numeric vector
#'   \item `vars_of_interest` must be a non-empty unique character vector
#'   \item `setequal(names(true_params), vars_of_interest)` must be TRUE
#'   \item No duplicate names in any named vector/list
#'   \item `meta` must be a named list with scalar values only (if present)
#' }
#'
#' @section Errors:
#' Throws a `bayesim_data_error` condition if validation fails.
#'
#' @keywords internal
#' @export
#'
#' @examples
#' # Valid data bundle
#' data_bundle <- list(
#'   train = data.frame(x = 1:10, y = rnorm(10)),
#'   test = data.frame(x = 11:15, y = rnorm(5)),
#'   response = "y",
#'   true_params = c(beta = 1.5, sigma = 0.5),
#'   vars_of_interest = c("beta", "sigma")
#' )
#' validate_data_bundle(data_bundle)
validate_data_bundle <- function(data_bundle) {
  if (!is.list(data_bundle)) {
    stop(
      bayesim_data_error(
        "data_bundle must be a list, got " %+% class(data_bundle)[1]
      )
    )
  }

  if (is.null(data_bundle$train)) {
    stop(bayesim_data_error("data_bundle$train is required and cannot be NULL"))
  }
  if (!is.data.frame(data_bundle$train) && !is.matrix(data_bundle$train)) {
    stop(bayesim_data_error(
      "data_bundle$train must be a data.frame or matrix, got " %+%
        class(data_bundle$train)[1]
    ))
  }
  if (nrow(data_bundle$train) < 1) {
    stop(bayesim_data_error("data_bundle$train must have at least 1 row"))
  }

  if (!is.null(data_bundle$test)) {
    if (
      !is.data.frame(data_bundle$test) &&
        !is.matrix(data_bundle$test) &&
        !is.list(data_bundle$test)
    ) {
      stop(bayesim_data_error(
        "data_bundle$test must be NULL or a structured object (data.frame, matrix, or list)"
      ))
    }
  }

  # response is required (a data bundle must name the outcome column).
  if (is.null(data_bundle$response)) {
    stop(bayesim_data_error(
      "data_bundle$response is required (the name of the response column)"
    ))
  }
  if (!is.character(data_bundle$response)) {
    stop(
      bayesim_data_error(
        "data_bundle$response must be a character vector, got " %+%
          typeof(data_bundle$response)
      )
    )
  }
  if (
    !is.null(data_bundle$response) &&
      (length(data_bundle$response) < 1 ||
        anyNA(data_bundle$response) ||
        any(data_bundle$response == ""))
  ) {
    stop(bayesim_data_error(
      "data_bundle$response cannot contain NA or empty strings"
    ))
  }

  if (
    !is.null(data_bundle$response) &&
      is.data.frame(data_bundle$train) &&
      length(data_bundle$response) == 1 &&
      !(data_bundle$response %in% names(data_bundle$train))
  ) {
    stop(
      bayesim_data_error(
        "data_bundle$response '" %+%
          data_bundle$response %+%
          "' not found in train columns: " %+%
          paste(names(data_bundle$train), collapse = ", ")
      )
    )
  }
  if (
    !is.null(data_bundle$test) &&
      !is.null(data_bundle$response) &&
      is.data.frame(data_bundle$test) &&
      length(data_bundle$response) == 1 &&
      !(data_bundle$response %in% names(data_bundle$test))
  ) {
    stop(
      bayesim_data_error(
        "data_bundle$response '" %+%
          data_bundle$response %+%
          "' not found in test columns: " %+%
          paste(names(data_bundle$test), collapse = ", ")
      )
    )
  }

  # E6: true_params and vars_of_interest are OPTIONAL (jointly NULL). Pure
  # model-comparison / predictive studies on truth-free data have no truths;
  # truth-dependent metrics already degrade to NA. Keep the integrity checks
  # when present.
  if (!is.null(data_bundle$true_params)) {
    if (!is.numeric(data_bundle$true_params)) {
      stop(
        bayesim_data_error(
          "data_bundle$true_params must be a numeric vector, got " %+%
            typeof(data_bundle$true_params)
        )
      )
    }
    if (is.null(names(data_bundle$true_params))) {
      stop(bayesim_data_error("data_bundle$true_params must have names"))
    }
    if (anyDuplicated(names(data_bundle$true_params)) > 0) {
      dup_names <- names(data_bundle$true_params)[duplicated(names(
        data_bundle$true_params
      ))]
      stop(
        bayesim_data_error(
          "data_bundle$true_params has duplicate names: " %+%
            paste(unique(dup_names), collapse = ", ")
        )
      )
    }
  }

  if (!is.null(data_bundle$vars_of_interest)) {
    if (!is.character(data_bundle$vars_of_interest)) {
      stop(
        bayesim_data_error(
          "data_bundle$vars_of_interest must be a character vector, got " %+%
            typeof(data_bundle$vars_of_interest)
        )
      )
    }
    if (length(data_bundle$vars_of_interest) < 1) {
      stop(bayesim_data_error("data_bundle$vars_of_interest must be non-empty"))
    }
    if (anyDuplicated(data_bundle$vars_of_interest) > 0) {
      dup_vars <- data_bundle$vars_of_interest[duplicated(
        data_bundle$vars_of_interest
      )]
      stop(
        bayesim_data_error(
          "data_bundle$vars_of_interest has duplicate values: " %+%
            paste(unique(dup_vars), collapse = ", ")
        )
      )
    }
    if (
      anyNA(data_bundle$vars_of_interest) ||
        any(data_bundle$vars_of_interest == "")
    ) {
      stop(bayesim_data_error(
        "data_bundle$vars_of_interest cannot contain NA or empty strings"
      ))
    }
  }

  # Validate true_params names match vars_of_interest (only when BOTH present).
  if (
    !is.null(data_bundle$true_params) &&
      !is.null(data_bundle$vars_of_interest) &&
      !setequal(names(data_bundle$true_params), data_bundle$vars_of_interest)
  ) {
    in_params_not_vars <- setdiff(
      names(data_bundle$true_params),
      data_bundle$vars_of_interest
    )
    in_vars_not_params <- setdiff(
      data_bundle$vars_of_interest,
      names(data_bundle$true_params)
    )
    msg <- "data_bundle: names(true_params) must exactly match vars_of_interest. "
    if (length(in_params_not_vars) > 0) {
      msg <- msg %+%
        "In true_params but not vars_of_interest: " %+%
        paste(in_params_not_vars, collapse = ", ") %+%
        ". "
    }
    if (length(in_vars_not_params) > 0) {
      msg <- msg %+%
        "In vars_of_interest but not true_params: " %+%
        paste(in_vars_not_params, collapse = ", ")
    }
    stop(bayesim_data_error(msg))
  }

  # Optional: validate meta is a named list of scalars
  if (!is.null(data_bundle$meta)) {
    if (!is.list(data_bundle$meta)) {
      stop(
        bayesim_data_error(
          "data_bundle$meta must be a list, got " %+% class(data_bundle$meta)[1]
        )
      )
    }
    if (length(data_bundle$meta) > 0) {
      meta_names <- names(data_bundle$meta)
      if (is.null(meta_names) || any(meta_names == "" | is.na(meta_names))) {
        stop(bayesim_data_error(
          "data_bundle$meta must have non-empty names for all elements"
        ))
      }
      if (anyDuplicated(meta_names) > 0) {
        dup_names <- meta_names[duplicated(meta_names)]
        stop(
          bayesim_data_error(
            "data_bundle$meta has duplicate names: " %+%
              paste(unique(dup_names), collapse = ", ")
          )
        )
      }
      # Check that all meta values are scalars
      for (nm in meta_names) {
        val <- data_bundle$meta[[nm]]
        if (!is.null(val) && length(val) != 1) {
          stop(
            bayesim_data_error(
              "data_bundle$meta$" %+%
                nm %+%
                " must be a scalar, got length " %+%
                length(val)
            )
          )
        }
      }
    }
  }

  invisible(data_bundle)
}

# =============================================================================
# Fit Result Validation
# =============================================================================

#' Validate Fit Result Interface
#'
#' Validates that a fit_result conforms to the bayesim_fit_result interface.
#' This is a wrapper around [validate_bayesim_fit_result()] that provides
#' clear error messages for contract violations.
#'
#' @param fit_result A bayesim_fit_result object to validate.
#'
#' @return The input `fit_result`, invisibly, if validation passes.
#'
#' @details
#' The fit_result must satisfy the bayesim_fit_result class requirements:
#' \itemize{
#'   \item Must inherit from "bayesim_fit_result" class
#'   \item If `success` is TRUE, `error` must be NULL
#'   \item If `success` is FALSE, `error` must be non-NULL
#'   \item `timing$total` must be a non-negative scalar numeric
#'   \item If `draws` is not NULL, it must be a matrix with column names
#'   \item `warnings` must be a character vector
#'   \item `diagnostics` must be a list
#' }
#'
#' @section Errors:
#' Throws a `bayesim_contract_error` condition if validation fails.
#'
#' @keywords internal
#'
#' @seealso [validate_bayesim_fit_result()], [new_fit_result()]
#'
#' @examples
#' \dontrun{
#' # Create a valid fit result
#' draws <- matrix(rnorm(100), ncol = 2, nrow = 50)
#' colnames(draws) <- c("alpha", "beta")
#' result <- new_fit_result(
#'   success = TRUE,
#'   draws = draws,
#'   diagnostics = list(rhat_max = 1.01),
#'   timing = list(total = 5.0, warmup = 2.5, sample = 2.5)
#' )
#' validate_fit_result_interface(result)
#' }
validate_fit_result_interface <- function(fit_result) {
  tryCatch(
    {
      validate_bayesim_fit_result(fit_result)
    },
    error = function(e) {
      stop(
        bayesim_contract_error(
          "Fit result contract violation: " %+% conditionMessage(e)
        )
      )
    }
  )
  invisible(fit_result)
}

# =============================================================================
# Fitter Interface Validation
# =============================================================================
# B3: the duplicate lightweight fitter class check (check_fitter_class /
# validate_fitter_interface) was merged into the exported validate_fitter()
# in R/fitter.R, which performs the class hierarchy + method checks. Internal
# callers use validate_fitter() directly.

# =============================================================================
# Metric Interface Validation
# =============================================================================

#' @title Validate a Metric Object
#'
#' @description
#' Validates that a metric object is an S7 instance of the Metric class with a
#' valid `name` property. This is the canonical metric validator (B3 merges the
#' former internal `validate_metric_interface` into this exported name). Method
#' existence is not checked because S7 dispatches via generics and the base
#' class raises errors for unimplemented abstract methods.
#'
#' @param metric An S7 object to validate as a Metric.
#'
#' @return The input `metric`, invisibly, if validation passes.
#'
#' @section Errors:
#' Throws a `bayesim_contract_error` condition if validation fails.
#'
#' @export
#'
#' @seealso [Metric], [validate_metric_output()], [validate_fitter()]
validate_metric <- function(metric) {
  if (!S7::S7_inherits(metric)) {
    stop(
      bayesim_contract_error(
        "metric must be an S7 object, got " %+% class(metric)[1]
      )
    )
  }

  if (!S7::S7_inherits(metric, Metric)) {
    stop(
      bayesim_contract_error(
        "metric must inherit from Metric class, got class: " %+%
          paste(class(metric), collapse = ", ")
      )
    )
  }

  metric_name <- tryCatch(
    {
      metric@name
    },
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(metric_name) ||
      !is.character(metric_name) ||
      length(metric_name) != 1
  ) {
    stop(
      bayesim_contract_error(
        "metric must have a 'name' property that is a scalar character"
      )
    )
  }

  if (is.na(metric_name) || metric_name == "") {
    stop(
      bayesim_contract_error(
        "metric name cannot be NA or empty string"
      )
    )
  }

  invisible(metric)
}

# =============================================================================
# Simulation Config Validation
# =============================================================================

#' Validate Simulation Configuration
#'
#' Validates that a SimulationConfig object is complete and properly configured
#' for running a simulation. This includes validating the fitter interface,
#' all metrics, and the data_generator function signature.
#'
#' @param config An S7 SimulationConfig object to validate.
#'
#' @return TRUE if validation passes.
#'
#' @details
#' The configuration is validated for:
#' \itemize{
#'   \item Being a valid SimulationConfig S7 object
#'   \item Having a non-NULL fitter that passes [validate_fitter()]
#'   \item Having metrics (if present) that pass [validate_metric()]
#'   \item Having a data_generator function with the correct signature
#' }
#'
#' @section Errors:
#' Throws a `bayesim_config_error` condition if validation fails.
#'
#' @keywords internal
#'
#' @seealso [simulation_config()], [validate_fitter()], [validate_metric()]
#'
#' @examples
#' \dontrun{
#' config <- simulation_config(
#'   data_grid = data.frame(n = 100),
#'   fit_grid = data.frame(model = "baseline"),
#'   data_generator = my_data_gen,
#'   fitter = my_fitter,
#'   metrics = list(my_metric),
#'   n_replicates = 10L,
#'   seed = 42L
#' )
#' validate_simulation_config(config)
#' }
validate_simulation_config <- function(config) {
  if (!is_simulation_config(config)) {
    stop(
      bayesim_config_error(
        "config must be a SimulationConfig object, got " %+% class(config)[1]
      )
    )
  }

  # Fitter is required for simulation
  if (is.null(config@fitter)) {
    stop(
      bayesim_config_error(
        "config@fitter cannot be NULL - a fitter is required to run simulations"
      )
    )
  }

  tryCatch(
    {
      validate_fitter(config@fitter)
    },
    error = function(e) {
      stop(
        bayesim_config_error(
          "Invalid fitter in configuration: " %+% conditionMessage(e)
        )
      )
    }
  )

  # Each metric must pass class/name validation
  if (!is.null(config@metrics) && length(config@metrics) > 0) {
    for (i in seq_along(config@metrics)) {
      metric <- config@metrics[[i]]
      tryCatch(
        {
          validate_metric(metric)
        },
        error = function(e) {
          metric_id <- if (S7::S7_inherits(metric) && !is.null(metric@name)) {
            metric@name
          } else if (is.list(metric) && !is.null(metric$name)) {
            metric$name
          } else {
            paste0("[", i, "]")
          }
          stop(
            bayesim_config_error(
              "Invalid metric " %+%
                metric_id %+%
                " in configuration: " %+%
                conditionMessage(e)
            )
          )
        }
      )
    }
  }

  # data_generator must be a function with at least 3 args
  if (is.null(config@data_generator)) {
    stop(
      bayesim_config_error(
        "config@data_generator cannot be NULL"
      )
    )
  }

  if (!is.function(config@data_generator)) {
    stop(
      bayesim_config_error(
        "config@data_generator must be a function"
      )
    )
  }

  gen_formals <- names(formals(config@data_generator))
  required_args <- c("data_spec", "task_ctx")

  if (length(gen_formals) < length(required_args)) {
    stop(
      bayesim_config_error(
        "data_generator must accept at least 2 arguments: (data_spec, task_ctx). " %+%
          "Got " %+%
          length(gen_formals) %+%
          " arguments: " %+%
          paste(gen_formals, collapse = ", ")
      )
    )
  }

  has_task_grid <- !is.null(config@task_grid) && nrow(config@task_grid) > 0
  has_crossed_grids <- !is.null(config@data_grid) &&
    nrow(config@data_grid) > 0 &&
    !is.null(config@fit_grid) &&
    nrow(config@fit_grid) > 0

  if (!has_task_grid && !has_crossed_grids) {
    stop(
      bayesim_config_error(
        "config must include task_grid or both data_grid and fit_grid"
      )
    )
  }

  if (
    is.null(config@n_replicates) ||
      length(config@n_replicates) != 1 ||
      is.na(config@n_replicates) ||
      config@n_replicates < 1
  ) {
    stop(
      bayesim_config_error(
        "config@n_replicates must be a positive integer >= 1"
      )
    )
  }

  if (is.null(config@seed) || length(config@seed) != 1 || is.na(config@seed)) {
    stop(
      bayesim_config_error(
        "config@seed must be a single non-NA integer"
      )
    )
  }

  TRUE
}

# =============================================================================
# Helper Functions
# =============================================================================
