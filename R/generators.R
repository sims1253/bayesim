# Generator library -------------------------------------------------------
#
# Factory constructors for simulation data generators:
#   - fixed_truth_generator():    user-supplied truth, fixed across tasks
#   - prior_predictive_generator(): draw theta from the model prior, then y|theta
#   - ifs_generator():            inverse forward sampling from a preconditioning fit
#
# All factories return a closure with the standard signature
#   (data_spec, task_ctx) -> data_bundle
# that consumes the AMBIENT RNG state (the worker restores the per-task
# L'Ecuyer stream via set_task_rng() before invoking the generator).
# `task_ctx$seed` carries an integer seed for backends (e.g. Stan) that need
# one, but generators must NOT use it to re-seed; doing so would defeat
# stream-based determinism.


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
#'   `test`, `meta`.
#'
#' @return A generator function `(data_spec, task_ctx) -> data_bundle`.
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

  function(data_spec, task_ctx) {
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
#' @return A generator function `(data_spec, task_ctx) -> data_bundle`.
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

  function(data_spec, task_ctx) {
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
#' @return A generator function `(data_spec, task_ctx) -> data_bundle`.
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

  function(data_spec, task_ctx) {
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

#' Default vars_of_interest for a brmsfit: population-level effects, plus the
#' residual scale `sigma` when present (E5: previously defaulted to population
#' effects only, silently excluding sigma/auxiliary parameters from SBC).
#' brms names effects "b_<coefname>"; strip the "b_" prefix for true_params names.
#' @keywords internal
.default_prior_vars <- function(fit) {
  vars <- brms::variables(fit)
  b_vars <- vars[grepl("^b_", vars)]
  out <- if (length(b_vars) > 0L) {
    sub("^b_", "", b_vars)
  } else {
    # Fall back to all variables if no population-level effects (rare).
    vars
  }
  # Include sigma (and other auxiliary scales) when present, so distributional
  # parameters are recovered/SBC-checked by default.
  aux <- vars[vars %in% c("sigma")]
  unique(c(out, aux))
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


# Fitter-agnostic generators --------------------------------------------------
#
# prior_draws_generator() and forward_sim_generator() are the fitter-agnostic
# counterparts of prior_predictive_generator() / ifs_generator(). Instead of
# requiring a brmsfit, they drive parameter/data generation through the S7
# Fitter interface (extract_draws(), predict_fit()), so they work with
# LinearRegressionFitter, BrmsFitter, CmdStanFitter, or any custom Fitter.
#
# Both factories fit ONCE on a pilot data_bundle supplied by the caller, store
# the resulting draws matrix, and then on each task pick a deterministic draw
# (indexed by task_ctx$rep_idx, wrapped modulo n_draws) and forward-simulate y
# for the predictors in data_spec.
#
# IMPORTANT (limitation): a true model prior is only directly available for
# brms (via brms::prior_draws). For other fitters there is no generic
# "prior_draws" S7 method, so prior_draws_generator() uses the DRAWS STORED ON
# THE FIT — for LinearRegressionFitter these are NIG-prior-posterior draws from
# the pilot (with a weak default prior, prior_precision = 1e-6), which is an
# APPROXIMATION of the prior predictive (concentrated on the pilot's posterior
# region, like IFS). For BrmsFitter, prior_draws_generator() first tries
# brms::prior_draws(); if unavailable, it degrades to the fit's posterior
# draws. Brms users who want the full prior predictive should use the
# brms-specific prior_predictive_generator().


#' Construct a fitter-agnostic prior-draws data generator
#'
#' `prior_draws_generator()` is the fitter-agnostic analogue of
#' [prior_predictive_generator()]. It works through the S7 [Fitter] interface
#' rather than brms-specific functions, so it can be used with
#' [LinearRegressionFitter], [BrmsFitter], [CmdStanFitter], or any custom
#' Fitter.
#'
#' The factory fits the model ONCE on `pilot_bundle` (provided by the caller),
#' extracts parameter draws via [extract_draws()], and stores them. The returned
#' closure, on each call, picks a draw deterministically indexed by
#' `task_ctx$rep_idx` (wrapped modulo the number of stored draws), uses it as
#' `true_params`, and forward-simulates `y` from the supplied predictors.
#'
#' Forward simulation: the response is drawn from a Gaussian with mean equal to
#' the linear predictor `X theta` (using the coefficient columns of the draw)
#' and standard deviation equal to the `sigma` column of the draw (if present,
#' else 1). This is the natural data-generating process for Gaussian linear
#' models — the common case for [LinearRegressionFitter].
#'
#' @section Limitations (non-brms):
#' A true model prior is directly accessible only for brms fits
#' ([brms::prior_draws()]). For other fitters there is no generic
#' "prior_draws" S7 method, so this factory falls back to the draws stored on
#' the pilot fit. For [LinearRegressionFitter] those are NIG-prior-conditioned
#' posterior draws (the prior is weak by default, `prior_precision = 1e-6`),
#' which makes this an *approximate* prior-predictive path concentrated on the
#' pilot's posterior region. Brms users who need full prior-predictive coverage
#' should prefer the brms-specific [prior_predictive_generator()].
#'
#' @param fitter An S7 [Fitter] object (e.g. [LinearRegressionFitter]).
#' @param fit_spec A list (single-row fit_grid entry) carrying at least
#'   `formula` (a base R formula). For [LinearRegressionFitter] a formula like
#'   `y ~ x` is expected.
#' @param pilot_bundle A `data_bundle` list (`train`, `response`, etc.) used for
#'   the one-time preconditioning fit. The caller is responsible for providing a
#'   representative pilot dataset. Must contain a `train` data.frame whose
#'   column names match the design implied by `fit_spec$formula`.
#' @param predictor_generator Function `(data_spec, task_ctx) -> data.frame`
#'   producing the design matrix of predictors (everything except the response).
#'   Must consume the ambient RNG state.
#' @param response Name of the response column. Defaults to
#'   `pilot_bundle$response`, falling back to the LHS of `fit_spec$formula`,
#'   then to `"y"`.
#' @param n_draws Optional integer override for the number of draws to store; if
#'   `NULL` (default), uses the number of draws returned by [extract_draws()].
#'
#' @return A generator function `(data_spec, task_ctx) -> data_bundle`.
#' @export
#' @seealso [prior_predictive_generator()], [forward_sim_generator()],
#'   [Fitter], [LinearRegressionFitter]
prior_draws_generator <- function(fitter, fit_spec, pilot_bundle,
                                  predictor_generator, response = NULL,
                                  n_draws = NULL) {
  if (!S7::S7_inherits(fitter, Fitter)) {
    stop(bayesim_config_error(
      "fitter must be an S7 Fitter object, got " %+% class(fitter)[1]
    ))
  }
  if (!is.function(predictor_generator)) {
    stop(bayesim_config_error("predictor_generator must be a function"))
  }
  if (is.null(pilot_bundle) || is.null(pilot_bundle$train)) {
    stop(bayesim_config_error(
      "pilot_bundle$train is required for the preconditioning fit"
    ))
  }

  resp <- response %||%
    pilot_bundle$response %||%
    .fit_spec_response_name(fit_spec) %||%
    "y"

  # BrmsFitter shortcut: if the fit is a sample_prior="only" brmsfit,
  # prefer true prior draws via brms::prior_draws() when available.
  draws_mat <- NULL
  prior_source <- "posterior_degraded"
  if (inherits(fitter, "BrmsFitter")) {
    fr <- .pilot_fit(fitter, pilot_bundle, fit_spec)
    prior_d <- tryCatch(
      {
        pd <- brms::prior_draws(fr$fit)
        if (is.null(pd) || nrow(pd) < 1L) NULL else pd
      },
      error = function(e) NULL
    )
    if (!is.null(prior_d)) {
      # Coerce prior_draws (a data.frame, one column per prior parameter) to a
      # draws matrix. Strip the "b_" prefix from population-level effects so
      # names line up with the cleaned convention used elsewhere.
      draws_mat <- as.matrix(prior_d)
      cn <- colnames(draws_mat)
      cleaned <- sub("^b_", "", cn)
      colnames(draws_mat) <- cleaned
      prior_source <- "prior"
    } else {
      draws_mat <- extract_draws(fitter, fr)
      prior_source <- "posterior_degraded"
    }
  } else {
    fr <- .pilot_fit(fitter, pilot_bundle, fit_spec)
    draws_mat <- extract_draws(fitter, fr)
  }

  if (is.null(draws_mat) || nrow(draws_mat) < 1L) {
    stop(bayesim_config_error(
      "pilot fit produced no usable draws for prior_draws_generator"
    ))
  }
  if (!is.null(n_draws)) {
    n_draws <- as.integer(n_draws)
    if (n_draws > 0L && n_draws < nrow(draws_mat)) {
      draws_mat <- draws_mat[seq_len(n_draws), , drop = FALSE]
    }
  }

  # vars_of_interest = the fitter's draw column names. Coefficient + scale
  # columns are reported as the truth; forward simulation uses them directly.
  voi <- colnames(draws_mat)
  if (is.null(voi) || length(voi) < 1L) {
    stop(bayesim_config_error(
      "pilot fit's draws matrix has no column names; cannot name true_params"
    ))
  }
  stored_n <- nrow(draws_mat)

  function(data_spec, task_ctx) {
    rep_idx <- task_ctx$rep_idx %||% 1L
    draw_id <- ((as.integer(rep_idx) - 1L) %% stored_n) + 1L

    newdata <- predictor_generator(data_spec, task_ctx)

    theta <- draws_mat[draw_id, , drop = FALSE]
    theta_vec <- as.numeric(theta)
    names(theta_vec) <- voi

    # Forward-simulate y from theta. Coefficient columns are everything except
    # sigma; the linear predictor is built by multiplying the design columns by
    # their coefficients (intercept column "Intercept" handled implicitly — it
    # is just another additive column with value 1 in the design matrix).
    sigma_col <- "sigma"
    coef_names <- setdiff(voi, sigma_col)
    mu <- .linear_predictor(newdata, theta_vec, coef_names, resp)
    sigma_val <- if (sigma_col %in% voi) theta_vec[[sigma_col]] else 1
    y <- mu + stats::rnorm(length(mu)) * sigma_val

    train <- newdata
    train[[resp]] <- y

    list(
      train = train,
      test = NULL,
      response = resp,
      true_params = theta_vec,
      vars_of_interest = voi,
      meta = list(
        generator = "prior_draws",
        truth_draw_id = draw_id,
        prior_source = prior_source
      )
    )
  }
}


#' Construct a fitter-agnostic forward-simulation (IFS) data generator
#'
#' `forward_sim_generator()` is the fitter-agnostic analogue of
#' [ifs_generator()]. It fits the model ONCE on `pilot_bundle`, draws theta from
#' the posterior (via [extract_draws()]), and forward-simulates `y` via the
#' fitter's [predict_fit()] (posterior-predictive) at the chosen draw's
#' parameters. Works with any Fitter that supports [predict_fit()]
#' (LinearRegressionFitter, BrmsFitter, ...).
#'
#' Unlike [prior_draws_generator()] (which targets the prior), this generator
#' concentrates the truth draw in a region of high posterior mass — the
#' canonical SBC generator for models with diffuse or improper priors.
#'
#' Because forward simulation here relies on [predict_fit()], the response is
#' drawn exactly as the fitter implements its posterior-predictive sampling,
#' which respects the fitter's response distribution and link function. Each
#' task uses a distinct draw, deterministically indexed by `task_ctx$rep_idx`.
#'
#' @section Limitations:
#' The `true_params` reported by this generator are the [extract_draws()]
#' columns for the selected draw; for [LinearRegressionFitter] these are
#' `Intercept`, `<coef>`, `sigma`. The response is forward-simulated via the
#' fitter's [predict_fit()] applied to the predictor design (a single Gaussian
#' draw for Gaussian fitters). Fitters without [predict_fit()] support are not
#' supported (e.g. raw [CmdStanFitter], which has no newdata semantics).
#'
#' @param fitter An S7 [Fitter] object supporting [predict_fit()].
#' @param fit_spec A list (single-row fit_grid entry) with at least `formula`.
#' @param pilot_bundle A `data_bundle` list used for the one-time
#'   preconditioning fit.
#' @param predictor_generator Function `(data_spec, task_ctx) -> data.frame`
#'   producing predictor covariates. Must consume the ambient RNG state.
#' @param response Name of the response column. Defaults to
#'   `pilot_bundle$response`, falling back to the LHS of `fit_spec$formula`,
#'   then to `"y"`.
#' @param n_draws Optional integer override for the number of draws to store.
#'
#' @return A generator function `(data_spec, task_ctx) -> data_bundle`.
#' @export
#' @seealso [ifs_generator()], [prior_draws_generator()], [Fitter],
#'   [LinearRegressionFitter]
forward_sim_generator <- function(fitter, fit_spec, pilot_bundle,
                                  predictor_generator, response = NULL,
                                  n_draws = NULL) {
  if (!S7::S7_inherits(fitter, Fitter)) {
    stop(bayesim_config_error(
      "fitter must be an S7 Fitter object, got " %+% class(fitter)[1]
    ))
  }
  if (!is.function(predictor_generator)) {
    stop(bayesim_config_error("predictor_generator must be a function"))
  }
  if (is.null(pilot_bundle) || is.null(pilot_bundle$train)) {
    stop(bayesim_config_error(
      "pilot_bundle$train is required for the preconditioning fit"
    ))
  }

  resp <- response %||%
    pilot_bundle$response %||%
    .fit_spec_response_name(fit_spec) %||%
    "y"

  fr <- .pilot_fit(fitter, pilot_bundle, fit_spec)
  draws_mat <- extract_draws(fitter, fr)
  if (is.null(draws_mat) || nrow(draws_mat) < 1L) {
    stop(bayesim_config_error(
      "pilot fit produced no usable draws for forward_sim_generator"
    ))
  }
  if (!is.null(n_draws)) {
    n_draws <- as.integer(n_draws)
    if (n_draws > 0L && n_draws < nrow(draws_mat)) {
      draws_mat <- draws_mat[seq_len(n_draws), , drop = FALSE]
    }
  }

  voi <- colnames(draws_mat)
  if (is.null(voi) || length(voi) < 1L) {
    stop(bayesim_config_error(
      "pilot fit's draws matrix has no column names; cannot name true_params"
    ))
  }
  stored_n <- nrow(draws_mat)

  function(data_spec, task_ctx) {
    rep_idx <- task_ctx$rep_idx %||% 1L
    draw_id <- ((as.integer(rep_idx) - 1L) %% stored_n) + 1L

    newdata <- predictor_generator(data_spec, task_ctx)

    theta_vec <- as.numeric(draws_mat[draw_id, , drop = FALSE])
    names(theta_vec) <- voi

    # Forward-simulate y from theta via the fitter's predict_fit(). We build a
    # single-draw fit_result slice so predict_fit() produces a draw for this
    # theta. For LinearRegressionFitter (and any fitter whose predict_fit reads
    # fit_result$draws) we construct the slice directly; for BrmsFitter we
    # forward to brms::posterior_predict with a single draw index.
    y <- .forward_sim_y(fitter, fr, newdata, draw_id, theta_vec, voi, resp)

    train <- newdata
    train[[resp]] <- y

    list(
      train = train,
      test = NULL,
      response = resp,
      true_params = theta_vec,
      vars_of_interest = voi,
      meta = list(
        generator = "forward_sim",
        truth_draw_id = draw_id
      )
    )
  }
}


# Internal helpers for the fitter-agnostic generators ----------------------

#' Run the one-time preconditioning pilot fit via fit_model().
#'
#' Uses the fitter's fit_model() with a fixed seed so the preconditioning fit is
#' reproducible across tasks (the worker's ambient RNG stream is independent of
#' the seed passed here, which only governs the MCMC/NIG draw generation for the
#' pilot). Returns the bayesim_fit_result.
#' @keywords internal
.pilot_fit <- function(fitter, pilot_bundle, fit_spec) {
  seed <- pilot_bundle$.pilot_seed %||% 0L
  fit_model(
    fitter,
    data_bundle = pilot_bundle,
    fit_spec = fit_spec,
    seed = seed,
    task_ctx = list(task_id = "pilot", rep_idx = 1L)
  )
}

#' Resolve the response name from a fit_spec (LHS of fit_spec$formula).
#' Returns NA_character_ if it cannot be resolved.
#' @keywords internal
.fit_spec_response_name <- function(fit_spec) {
  out <- tryCatch(
    {
      fml <- fit_spec$formula
      if (is.null(fml)) NA_character_ else all.vars(fml)[1L]
    },
    error = function(e) NA_character_
  )
  if (is.na(out) || !nzchar(out)) NA_character_ else out
}

#' Compute the linear predictor X theta for a design data.frame and a draw.
#'
#' coef_names names the additive coefficient columns; the intercept is named
#' "Intercept" and contributes a constant 1 (since predictor designs do not
#' include it as a column). Remaining coef_names must match columns of newdata.
#' Returns a length-nrow(newdata) numeric vector.
#' @keywords internal
.linear_predictor <- function(newdata, theta_vec, coef_names, resp) {
  n <- nrow(newdata)
  if (is.null(n) || n < 1L) {
    stop(bayesim_config_error(
      "predictor_generator returned an empty or invalid data.frame"
    ))
  }
  mu <- rep(0, n)
  for (cn in coef_names) {
    val <- theta_vec[[cn]]
    if (cn == "Intercept") {
      mu <- mu + val
    } else {
      if (!cn %in% names(newdata)) {
        stop(bayesim_config_error(
          "predictor_generator is missing column '" %+% cn %+%
            "' required by the stored draw"
        ))
      }
      mu <- mu + val * newdata[[cn]]
    }
  }
  mu
}

#' Forward-simulate the response for a single draw across fitters.
#'
#' Dispatches by fitter class:
#'   - LinearRegressionFitter (and fitters that expose fit_result$draws):
#'     reconstruct a one-draw fit_result and call predict_fit(); since
#'     predict_fit consumes the ambient RNG state for noise, this is a single
#'     fresh Gaussian draw at the selected theta.
#'   - BrmsFitter: call brms::posterior_predict(fit, newdata, draw_ids) for a
#'     single draw.
#'   - Fallback: try predict_fit() with a sliced fit_result; error on failure.
#'
#' Returns a numeric vector of length nrow(newdata).
#' @keywords internal
.forward_sim_y <- function(fitter, fit_result, newdata, draw_id, theta_vec,
                           voi, resp) {
  n <- nrow(newdata)

  if (inherits(fitter, "BrmsFitter")) {
    pp <- tryCatch(
      as.numeric(brms::posterior_predict(
        fit_result$fit,
        newdata = newdata,
        draw_ids = draw_id
      )),
      error = function(e) NULL
    )
    if (is.null(pp) || length(pp) != n) {
      stop(bayesim_config_error(
        "BrmsFitter posterior_predict failed during forward_sim_generator"
      ))
    }
    return(pp)
  }

  # Generic path: slice fit_result to the selected draw and ask predict_fit().
  # predict_fit's seed argument is left NULL so it uses the ambient RNG state.
  # Some fitters (e.g. LinearRegressionFitter) build their design matrix from
  # newdata + formula via model.frame(), which requires the response column to
  # be present even though predict_fit only consumes X. Inject a zero
  # placeholder for resp when missing; predict_fit ignores the response values.
  pp_newdata <- newdata
  if (!is.null(resp) && !(resp %in% names(pp_newdata))) {
    pp_newdata[[resp]] <- 0
  }
  sliced <- fit_result
  if (!is.null(fit_result$draws)) {
    sliced$draws <- fit_result$draws[draw_id, , drop = FALSE]
  }
  preds <- predict_fit(fitter, sliced, newdata = pp_newdata, seed = NULL)
  if (is.null(preds) || is.null(preds$predicted_samples)) {
    stop(bayesim_config_error(
      "fitter's predict_fit() returned NULL; fitter does not support " %+%
        "posterior-predictive sampling for forward_sim_generator"
    ))
  }
  ps <- preds$predicted_samples
  # predict_fit returns S x N (draws x observations); we sliced to one draw so
  # S == 1. Take the first (and only) row.
  if (is.matrix(ps)) {
    if (nrow(ps) < 1L || ncol(ps) != n) {
      stop(bayesim_config_error(
        "predict_fit() returned a predicted_samples matrix with wrong " %+%
          "dimensions during forward_sim_generator"
      ))
    }
    return(as.numeric(ps[1L, ]))
  }
  as.numeric(ps)
}
