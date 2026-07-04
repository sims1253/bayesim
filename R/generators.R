# Generator library -------------------------------------------------------
#
# Factory constructors for simulation data generators:
#   - fixed_truth_generator():    user-supplied truth, fixed across tasks
#   - prior_predictive_generator(): draw theta from the model prior, then y|theta
#   - ifs_generator():            inverse forward sampling from a preconditioning fit
#
# All factories return a closure with the standard signature
#   (data_spec, seed, task_ctx) -> data_bundle
# that consumes the AMBIENT RNG state (the worker restores the per-task
# L'Ecuyer stream via set_task_rng() before invoking the generator). The `seed`
# argument is retained for backends requiring an explicit integer but is NOT
# used to re-seed; doing so would defeat stream-based determinism.


# Fixed-truth generator ---------------------------------------------------

#' Construct a fixed-truth data generator
#'
#' Returns a generator closure that delegates to a user-supplied `draw_data`
#' function, pinning `data_bundle$true_params` to the supplied `truth` for every
#' task. Use this when the data-generating truth is known and constant (e.g. a
#' fixed effect size sweep).
#'
#' The `draw_data` function must have signature `(data_spec, task_ctx) ->
#' list(train, test, response, ...)` WITHOUT `true_params`/`vars_of_interest` —
#' those are injected by the factory from `truth`. It must consume the ambient
#' RNG state (do not call `set.seed`/`withr::with_seed` inside).
#'
#' @param truth Named numeric vector; the data-generating parameters, used as
#'   `true_params` and `vars_of_interest` for every generated bundle.
#' @param draw_data Function with signature `(data_spec, task_ctx)` returning a
#'   list with at least `train` (data.frame). May optionally return `response`,
#'   `test`, `references`, `meta`.
#'
#' @return A generator function `(data_spec, seed, task_ctx) -> data_bundle`.
#' @export
#' @examples
#' \dontrun{
#' gen <- fixed_truth_generator(
#'   truth = c(beta = 1, sigma = 1),
#'   draw_data = function(data_spec, task_ctx) {
#'     n <- data_spec$n %||% 20L
#'     x <- stats::rnorm(n)
#'     y <- x + stats::rnorm(n)
#'     list(train = data.frame(y = y, x = x), response = "y")
#'   }
#' )
#' }
fixed_truth_generator <- function(truth, draw_data) {
  if (!is.numeric(truth) || is.null(names(truth))) {
    stop(bayesim_config_error(
      "truth must be a named numeric vector, got "
        %+% class(truth)[1]
        %+% " with "
        %+% if (is.null(names(truth))) "no names" else "names"
    ))
  }
  if (!is.function(draw_data)) {
    stop(bayesim_config_error("draw_data must be a function"))
  }
  # Loosely check the signature; allow trailing ...
  fmls <- names(formals(draw_data))
  if (!all(c("data_spec", "task_ctx") %in% fmls)) {
    stop(bayesim_config_error(
      "draw_data must have arguments data_spec and task_ctx"
    ))
  }

  function(data_spec, seed, task_ctx) {
    bundle <- draw_data(data_spec, task_ctx)
    bundle$true_params <- truth
    bundle$vars_of_interest <- names(truth)
    if (is.null(bundle$meta)) bundle$meta <- list()
    bundle$meta$generator <- "fixed_truth"
    bundle
  }
}


# Prior-predictive generator ----------------------------------------------

#' Construct a prior-predictive data generator
#'
#' Draws a parameter vector theta from the model prior (via a `sample_prior =
#' "only"` brmsfit) and simulates data `y ~ p(y | theta)` using
#' [brms::posterior_predict()]. Each task uses a distinct prior draw,
#' deterministically indexed by `task_ctx$rep_idx`, so Simulation-Based
#' Calibration ranks are well-defined and resume is reproducible.
#'
#' The prior model is compiled once (a `sample_prior = "only"` brmsfit); reuse
#' it across tasks via the model bank or by passing the same object. Predictor
#' covariates not implied by the prior are supplied by `predictor_generator`.
#'
#' @param prior_fit A brmsfit compiled with `sample_prior = "only"` (or a
#'   formula + family + prior combination to be compiled; see Details). Must
#'   contain prior draws.
#' @param predictor_generator Function `(data_spec, task_ctx) -> data.frame`
#'   producing the design matrix of predictors (everything except the response).
#'   Must consume the ambient RNG state. If `NULL`, the prior_fit's own data is
#'   reused at its original size.
#' @param vars_of_interest Character vector naming the prior parameters to
#'   report as `true_params` (defaults to all population-level effects
#'   `"b_<name>"`, renamed to `<name>`).
#' @param response Name of the response column (defaults to the LHS of
#'   `prior_fit`'s formula).
#'
#' @return A generator function `(data_spec, seed, task_ctx) -> data_bundle`.
#' @export
prior_predictive_generator <- function(prior_fit,
                                       predictor_generator = NULL,
                                       vars_of_interest = NULL,
                                       response = NULL) {
  if (!inherits(prior_fit, "brmsfit")) {
    stop(bayesim_config_error(
      "prior_fit must be a brmsfit compiled with sample_prior = 'only'"
    ))
  }
  n_draws <- posterior::ndraws(prior_fit)
  if (n_draws < 1L) {
    stop(bayesim_config_error(
      "prior_fit contains no draws; compile with sample_prior = 'only'"
    ))
  }

  resp <- response %||% .fit_response_name(prior_fit)
  voi <- vars_of_interest %||% .default_prior_vars(prior_fit)

  function(data_spec, seed, task_ctx) {
    # Deterministic draw index from the replicate index, wrapped into range.
    # rep_idx is 1-based; cycle modulo the number of prior draws so a study
    # with more replicates than prior draws still covers them all.
    rep_idx <- task_ctx$rep_idx %||% 1L
    draw_id <- ((as.integer(rep_idx) - 1L) %% n_draws) + 1L

    newdata <- if (is.function(predictor_generator)) {
      predictor_generator(data_spec, task_ctx)
    } else {
      prior_fit$data
    }

    # Forward-sample the response from the selected prior draw.
    y <- as.numeric(brms::posterior_predict(
      prior_fit,
      newdata = newdata,
      draw_ids = draw_id
    ))

    # Extract theta = the drawn parameters (population-level effects + sigma).
    draws_mat <- posterior::as_draws_matrix(prior_fit)
    theta <- .extract_truth(draws_mat, draw_id, voi)

    train <- newdata
    train[[resp]] <- y

    list(
      train = train,
      test = NULL,
      response = resp,
      true_params = theta,
      vars_of_interest = names(theta),
      references = setNames(rep(NA_real_, length(theta)), names(theta)),
      meta = list(generator = "prior_predictive", truth_draw_id = draw_id)
    )
  }
}


# IFS generator -----------------------------------------------------------

#' Construct an inverse forward sampling (IFS) data generator
#'
#' Draws a parameter vector theta from a preconditioning fit (typically a
#' posterior fit on a pilot dataset), uses it as the data-generating truth, and
#' forward-simulates data `y ~ p(y | theta)` respecting multivariate response
#' dependencies. Each task uses a distinct draw, deterministically indexed by
#' `task_ctx$rep_idx`.
#'
#' Unlike [prior_predictive_generator()] (which samples from the prior), IFS
#' samples from a preconditioning posterior, concentrating the truth draw in a
#' region of high posterior mass. This is the canonical SBC generator for
#' models with diffuse or improper priors.
#'
#' @param prefit A brmsfit with posterior draws to sample theta from (the
#'   preconditioning fit).
#' @param predictor_generator Function `(data_spec, task_ctx) -> data.frame`
#'   producing predictor covariates. Must consume the ambient RNG state. If
#'   `NULL`, `prefit$data` is reused.
#' @param vars_of_interest Character vector naming the parameters to report as
#'   `true_params`. Defaults to population-level effects.
#' @param response Name of the response column; defaults to the LHS of
#'   `prefit`'s formula.
#' @param lower_bound,upper_bound Optional numeric bounds for the response
#'   domain. If both `NULL` (default), no bounds are applied. If set, the
#'   out-of-bounds policy is governed by `truncate`:
#'   - `truncate = FALSE` (default): out-of-bounds draws set the response to
#'     all `NA`, which fails downstream validation and drops the replicate.
#'     NOTE: this is NOT a draw-level resample — it drops the replicate's draw
#'     from the SBC rank distribution, which biases ranks when violations are
#'     non-uniform across the posterior. Use only when the bounds are soft.
#'   - `truncate = TRUE`: clamp out-of-bounds response values to the nearest
#'     bound. This also biases the rank distribution (toward the bounds) but
#'     keeps the replicate.
#'   Neither option implements the plan's deterministic draw-resampling; both
#'   are documented honestly here. A future version may add
#'   `on_violation = "resample"`.
#' @param truncate Logical; if `TRUE` and bounds are set, clamp out-of-bounds
#'   response values instead of producing NAs. Default `FALSE`.
#'
#' @return A generator function `(data_spec, seed, task_ctx) -> data_bundle`.
#' @export
ifs_generator <- function(prefit,
                          predictor_generator = NULL,
                          vars_of_interest = NULL,
                          response = NULL,
                          lower_bound = NULL,
                          upper_bound = NULL,
                          truncate = FALSE) {
  if (!inherits(prefit, "brmsfit")) {
    stop(bayesim_config_error("prefit must be a brmsfit object"))
  }
  n_draws <- posterior::ndraws(prefit)
  if (n_draws < 1L) {
    stop(bayesim_config_error(
      "prefit contains no posterior draws; fit with at least one chain"
    ))
  }

  use_bounds <- !is.null(lower_bound) || !is.null(upper_bound)
  resp <- response %||% .fit_response_name(prefit)
  voi <- vars_of_interest %||% .default_prior_vars(prefit)

  function(data_spec, seed, task_ctx) {
    rep_idx <- task_ctx$rep_idx %||% 1L
    draw_id <- ((as.integer(rep_idx) - 1L) %% n_draws) + 1L

    newdata <- if (is.function(predictor_generator)) {
      predictor_generator(data_spec, task_ctx)
    } else {
      prefit$data
    }

    # Forward-sample the response(s) from the selected draw, respecting
    # multivariate dependency order. For univariate models this is a single
    # posterior_predict call; for multivariate brms models the topological
    # ordering in brms_response_sequence() ensures dependent responses see
    # their already-simulated predecessors. brms_full_ppred() now takes a
    # SINGLE draw and returns a single data.frame, eliminating the previous
    # draw-value-as-list-index mismatch (F1 root cause #3).
    sim_frame <- brms_full_ppred(prefit, newdata = newdata, draw = draw_id)

    # Bounds policy.
    if (use_bounds && resp %in% names(sim_frame)) {
      if (isTRUE(truncate)) {
        if (!is.null(lower_bound)) {
          sim_frame[[resp]] <- pmax(sim_frame[[resp]], lower_bound)
        }
        if (!is.null(upper_bound)) {
          sim_frame[[resp]] <- pmin(sim_frame[[resp]], upper_bound)
        }
      } else {
        y <- sim_frame[[resp]]
        oob <- FALSE
        if (!is.null(lower_bound) && any(y < lower_bound, na.rm = TRUE)) oob <- TRUE
        if (!is.null(upper_bound) && any(y > upper_bound, na.rm = TRUE)) oob <- TRUE
        if (oob) {
          # Out-of-bounds policy: NA-out the response (NOT a draw-level
          # resample). This drops the replicate, biasing SBC ranks when
          # violations are non-uniform — see ?ifs_generator. Deterministic
          # because draw_id is fixed by rep_idx.
          sim_frame[[resp]] <- NA_real_
        }
      }
    }

    draws_mat <- posterior::as_draws_matrix(prefit)
    theta <- .extract_truth(draws_mat, draw_id, voi)

    list(
      train = sim_frame,
      test = NULL,
      response = resp,
      true_params = theta,
      vars_of_interest = names(theta),
      references = setNames(rep(NA_real_, length(theta)), names(theta)),
      meta = list(generator = "ifs", truth_draw_id = draw_id)
    )
  }
}


# Internal helpers --------------------------------------------------------

#' Extract the response variable name from a brmsfit.
#' Falls back through formula forms since sample_prior="only" fits may strip
#' the standard formula() variables.
#' @keywords internal
.fit_response_name <- function(fit) {
  # Try the stored brmsformula first (most reliable).
  resp <- tryCatch(
    all.vars(fit$formula$formula)[1],
    error = function(e) NA_character_
  )
  if (!is.na(resp) && nchar(resp) > 0L) return(resp)
  # Fall back to the names of the original data.
  if (!is.null(fit$data) && ncol(fit$data) > 0L) {
    return(names(fit$data)[1L])
  }
  "y"
}

#' Default vars_of_interest for a brmsfit: population-level effects.
#' brms names them "b_<coefname>"; strip the "b_" prefix for true_params names.
#' @keywords internal
.default_prior_vars <- function(fit) {
  vars <- brms::variables(fit)
  b_vars <- vars[grepl("^b_", vars)]
  if (length(b_vars) == 0L) {
    # Fall back to all variables if no population-level effects (rare).
    return(vars)
  }
  # Strip the "b_" prefix, but keep the intercept distinguishable.
  out <- sub("^b_", "", b_vars)
  # Rename the intercept: brms stores it as "b_Intercept" (dense) or
  # "b_<resp>_Intercept" (sparse); keep the cleaned form.
  out
}

#' Extract a named parameter vector from a draws matrix at a given draw index.
#' Errors if a requested variable cannot be found (neither as the cleaned name
#' nor as "b_<name>"), since a silent NA would corrupt downstream SBC ranks.
#' @keywords internal
.extract_truth <- function(draws_mat, draw_id, vars_of_interest) {
  available <- colnames(draws_mat)
  candidates <- paste0("b_", vars_of_interest)
  hits <- vars_of_interest
  names(hits) <- vars_of_interest
  missing <- character(0)
  for (i in seq_along(vars_of_interest)) {
    v <- vars_of_interest[i]
    if (v %in% available) {
      next
    } else if (candidates[i] %in% available) {
      hits[i] <- candidates[i]
    } else {
      missing <- c(missing, v)
    }
  }
  if (length(missing) > 0L) {
    stop(bayesim_config_error(
      "vars_of_interest not found in model draws: "
        %+% paste(missing, collapse = ", ")
        %+% ". Available: "
        %+% paste(utils::head(available, 20), collapse = ", ")
    ))
  }
  vals <- draws_mat[draw_id, hits, drop = FALSE]
  out <- as.numeric(vals)
  names(out) <- vars_of_interest
  out
}


# Forward-sampling internals (ported from 0.x, decoupled from SBC/bayeshear)
# -------------------------------------------------------------------------

#' Full forward sampling of a brms model's responses for a SINGLE posterior draw
#'
#' Simulates the response(s) of a brmsfit at a single posterior draw, respecting
#' dependency order among multiple responses. Ported from the SBC package and
#' the 0.x bayesim codebase; bayeshear/SBC/future dependencies removed.
#'
#' `ifs_generator()` only ever needs one draw per task (the deterministic
#' `rep_idx`-indexed draw), so this function takes a single `draw` index and
#' returns a single data.frame. This eliminates the prior index mismatch where
#' results were stored at `pp_data[[draw_value]]` but read at `simulated[[1]]`.
#'
#' @param fit A brmsfit with posterior draws.
#' @param newdata Optional data.frame of predictors. If `NULL`, uses `fit$data`.
#' @param draw Single integer draw index to simulate.
#'
#' @return A data.frame (a copy of `newdata`) with the simulated response
#'   column(s) filled in.
#'
#' @keywords internal
brms_full_ppred <- function(fit, newdata = NULL, draw = 1L) {
  resp <- brms_response_sequence(fit)
  if (is.null(newdata)) {
    newdata <- fit$data
  }
  n <- nrow(newdata)
  pp_data <- newdata
  # Do NOT call brms::validate_newdata here: the response column is absent (we
  # are about to simulate it), and validate_newdata requires it. posterior_predict
  # performs its own internal validation per call.
  for (vars in resp) {
    pp_data[, vars] <- array(
      brms::posterior_predict(
        fit,
        newdata = pp_data,
        resp = vars,
        draw_ids = draw,
        allow_new_levels = TRUE
      ),
      dim = c(1, n, length(vars))
    )[1, , ]
  }
  pp_data
}

#' Topological sort of response variables by dependency depth (Kahn's
#' algorithm). Nodes with no incoming edges (predictors that depend on no other
#' response) come first.
#' @keywords internal
nodes_by_depth <- function(adj_matrix) {
  depth_list <- list()
  var_names <- rownames(adj_matrix)
  while (nrow(adj_matrix) > 0L) {
    pos <- which(apply(adj_matrix, 1, sum) == 0)
    if (length(pos) == 0L) {
      stop(bayesim_config_error(
        "Unsolvable variable dependencies detected: cyclic response structure. "
          %+% "The brms model has responses that depend on each other circularly; "
          %+% "forward sampling cannot order them."
      ))
    }
    depth_list <- c(depth_list, list(var_names[pos]))
    var_names <- var_names[-pos]
    adj_matrix <- adj_matrix[-pos, -pos, drop = FALSE]
  }
  depth_list
}

#' Determine the response dependency sequence of a brms model.
#'
#' S3 generic dispatched on the model object type. Returns a list of character
#' vectors: each element is the set of response names at a given dependency
#' depth, to be simulated together (topologically ordered so that a response
#' depending on another response is simulated after its predecessor). This is
#' the value iterated over by `brms_full_ppred()`'s `for (vars in resp)` loop,
#' so each element must be a *response-name* group, not a dependency list.
#'
#' Ported from the 0.x `brms_response_sequence.bform` method. Methods are
#' provided for both `brmsformula` (univariate) and `mvbrmsformula`
#' (multivariate); both feed into the adjacency-matrix + `nodes_by_depth()`
#' construction (the 0.x `bform` logic), so dispatch works regardless of
#' whether the `"bform"` class alias is registered.
#'
#' @keywords internal
brms_response_sequence <- function(x) {
  UseMethod("brms_response_sequence")
}

#' @method brms_response_sequence brmsfit
#' @exportS3Method
#' @keywords internal
brms_response_sequence.brmsfit <- function(x) {
  brms_response_sequence(x$formula)
}

#' @method brms_response_sequence brmsformula
#' @exportS3Method
#' @keywords internal
brms_response_sequence.brmsformula <- function(x) {
  # brmsterms() -> dependency term_list (named by response); then build the
  # adjacency matrix and return depth-ordered response-name groups.
  term_list <- brms_response_sequence(brms::brmsterms(x))
  .bform_response_depth(term_list)
}

#' @method brms_response_sequence mvbrmsformula
#' @exportS3Method
#' @keywords internal
brms_response_sequence.mvbrmsformula <- function(x) {
  # mvbrmsterms -> flattened dependency term_list across all responses; then
  # build the combined adjacency matrix and return depth-ordered groups.
  term_list <- brms_response_sequence(brms::brmsterms(x))
  .bform_response_depth(term_list)
}

# Shared core of the 0.x `brms_response_sequence.bform` method: given a named
# `term_list` (names = response names, values = character vector of variables
# each response depends on), build the response-dependency adjacency matrix and
# return `nodes_by_depth(adjacency)` — a list of response-name character
# vectors ordered by dependency depth.
.bform_response_depth <- function(term_list) {
  resp_vars <- names(term_list)
  adjacency <- t(sapply(term_list, function(x) is.element(resp_vars, x)))
  attr(adjacency, "dimnames") <- list(resp_vars, resp_vars)
  nodes_by_depth(adjacency)
}

#' @method brms_response_sequence mvbrmsterms
#' @exportS3Method
#' @keywords internal
brms_response_sequence.mvbrmsterms <- function(x) {
  # brms_response_sequence() over each brmsterms sub-object returns a named
  # list (names = response, values = dependency vars); sapply flattens the
  # list-of-length-1-lists into a single named term_list spanning all
  # responses, exactly as the 0.x bform method expected.
  names(x$terms) <- NULL
  term_list <- sapply(x$terms, brms_response_sequence)
  term_list
}

#' @method brms_response_sequence brmsterms
#' @exportS3Method
#' @keywords internal
brms_response_sequence.brmsterms <- function(x) {
  # Build the dependency term_list for a single response: a named list with one
  # element whose name is the response and whose value is the character vector
  # of variables the response depends on (via its dpars' btl formulas). For
  # bf(y ~ x) this is list(y = c("1", "x")). The formula-level method then
  # turns this into depth-ordered response-name groups.
  term_list <- list(unique(unlist(lapply(x$dpars, brms_response_sequence))))
  names(term_list) <- all.vars(x$respform)
  term_list
}

#' @method brms_response_sequence btl
#' @exportS3Method
#' @keywords internal
brms_response_sequence.btl <- function(x) {
  # 0.x semantics: the predictors this response depends on, prefixed with "1"
  # (intercept marker). For bf(y ~ x), all.vars(~x) = "x", so this returns
  # c("1", "x"); the adjacency row for "y" is then is.element(resp_vars,
  # c("1","x")) = all FALSE, giving "y" zero indegree (simulated first).
  c("1", all.vars(x$formula))
}

#' @method brms_response_sequence default
#' @exportS3Method
#' @keywords internal
brms_response_sequence.default <- function(x) {
  # Non-btl dpars (e.g. non-linear "btn" terms, or dpar-specific shape terms)
  # contribute no cross-response dependency. Returning character(0) lets the
  # dpars walk in brms_response_sequence.brmsterms skip them gracefully instead
  # of throwing an opaque "no applicable method" error. Inherited latent
  # fragility from the 0.x port; hardened here.
  character(0)
}
