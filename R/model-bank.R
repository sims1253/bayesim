#' @title Model Bank for BrmsFitter
#' @description Internal infrastructure that compiles each distinct brms model
#'   spec once (via `brms::brm(chains = 0)`) and reuses the resulting prefit
#'   across simulation tasks via `stats::update(recompile = FALSE)`. This
#'   eliminates catastrophic recompilation in studies that fit the same model
#'   to thousands of datasets.
#'
#'   The bank is transported to mirai daemons via a session option
#'   (`bayesim.model_bank`) rather than as an S7 fitter property, so it does not
#'   corrupt the config fingerprint (`capture_fitter_spec()` hashes S7
#'   properties only). See [set_model_bank()] / [get_model_bank()].
#'
#' @name model-bank
#' @keywords internal
NULL

# ============================================================================
# Session-level storage (mirai-serializable option)
# ============================================================================

#' Set the session model bank
#'
#' Stores the model bank in `options(bayesim.model_bank)`. Called on the
#' controller after the bank is built, and pushed to each daemon via
#' `mirai::everywhere()`. Pass NULL to clear it.
#'
#' @param bank A named list of `brmsfit` prefit objects keyed by
#'   [model_spec_hash()], or NULL.
#'
#' @return Invisible NULL. Called for its side effect.
#'
#' @keywords internal
set_model_bank <- function(bank) {
  options(bayesim.model_bank = bank)
  invisible(NULL)
}

#' Get the session model bank
#'
#' Retrieves the model bank set by [set_model_bank()]. Returns NULL when no bank
#' is active (e.g. `precompile = FALSE` or before [run_simulation()] builds it).
#'
#' @return A named list of `brmsfit` prefit objects, or NULL.
#'
#' @keywords internal
get_model_bank <- function() {
  getOption("bayesim.model_bank")
}

# ============================================================================
# Spec resolution and hashing
# ============================================================================

#' Resolve a family specification to a canonical brms family object
#'
#' Accepts a brms family object, a base `family` object, a string (resolved via
#' `brms::brmsfamily()`), or NULL (defaults to `brms::gaussian()`). Base
#' `family` objects (e.g. from `stats::gaussian()`) are canonicalized through
#' `brms::brmsfamily()` so that `gaussian()` and `brmsfamily("gaussian")` produce
#' identical objects (and thus identical model-spec hashes and Stan code).
#'
#' @param family A brms/base family object, a string, or NULL.
#'
#' @return A brms family object.
#'
#' @keywords internal
resolve_family <- function(family) {
  if (is.null(family)) {
    return(brms::brmsfamily("gaussian"))
  }
  if (is.character(family) && length(family) == 1L) {
    return(brms::brmsfamily(family))
  }
  # Canonicalize base family objects (stats::gaussian(), etc.) into brmsfamily
  # so the family identity (family+link) hashes identically to the string form.
  if (inherits(family, "family") && !inherits(family, "brmsfamily")) {
    fam_name <- family$family
    link_name <- family$link
    resolved <- tryCatch(
      brms::brmsfamily(fam_name, link = link_name),
      error = function(e) NULL
    )
    if (!is.null(resolved)) {
      return(resolved)
    }
  }
  family
}

#' Coerce a formula to a brmsformula
#'
#' Plain formulas are wrapped via `brms::::brmsformula()` equivalent
#' (`brms::bf()`); existing brmsformulas are returned unchanged.
#'
#' @param formula A formula or brmsformula.
#'
#' @return A brmsformula object.
#'
#' @keywords internal
resolve_formula <- function(formula) {
  if (inherits(formula, "brmsformula")) {
    return(formula)
  }
  brms::brmsformula(formula)
}

#' Compute a stable hash for a model spec
#'
#' Produces a stable string key identifying a brms model spec up to the fields
#' that affect Stan code generation: formula, family, prior, stanvars, and
#' backend. Uses `digest::digest()` on a normalized list representation, which
#' is deterministic across R sessions (unlike environment-bound object hashing).
#'
#' @param formula A formula, brmsformula, or string (resolved to brmsformula).
#' @param family A brms family object, string, or NULL (resolved).
#' @param prior A brms prior object or NULL.
#' @param stanvars A brms stanvars object or NULL.
#' @param backend Character scalar, e.g. "cmdstanr".
#'
#' @return A length-1 character string hash.
#'
#' @keywords internal
model_spec_hash <- function(formula, family, prior, stanvars, backend) {
  formula <- resolve_formula(formula)
  family <- resolve_family(family)

  # Normalize each component to a serializable representation. brms objects
  # carry environments that make rlang::hash() unstable across sessions, so we
  # digest a structured list instead.
  key <- list(
    formula = deparse_formula_recursive(formula),
    family = normalize_family(family),
    prior = normalize_prior(prior),
    stanvars = normalize_stanvars(stanvars),
    backend = backend
  )

  digest::digest(key, algo = "xxhash64")
}

#' Recursively deparse a brmsformula to a stable character representation
#'
#' @keywords internal
deparse_formula_recursive <- function(formula) {
  formula <- resolve_formula(formula)
  parts <- list(
    formula = paste(deparse(formula$formula), collapse = " "),
    pforms = if (!is.null(formula$pforms)) {
      paste(capture.output(print(formula$pforms)), collapse = " ")
    } else {
      ""
    },
    addition = if (!is.null(formula$addition)) {
      paste(capture.output(print(formula$addition)), collapse = " ")
    } else {
      ""
    }
  )
  parts
}

#' Normalize a brms family to a hashable list
#'
#' @keywords internal
normalize_family <- function(family) {
  family <- resolve_family(family)
  out <- list(
    family = family$family %||% NA_character_,
    link = family$link %||% NA_character_,
    class = class(family)[1]
  )
  # Include any extra brmsfamily flags (e.g. thresholds) if present as
  # non-environment slots. Bound vars/environments are intentionally dropped.
  extra <- family[
    names(family) %in%
      c(
        "threshold",
        "threshold_fixed",
        "nlpar",
        "dyad",
        "resp",
        "dpar"
      )
  ]
  if (length(extra) > 0L) {
    out$extra <- capture.output(print(extra))
  }
  out
}

#' Normalize a brms prior object to a hashable data.frame representation
#'
#' @keywords internal
normalize_prior <- function(prior) {
  if (is.null(prior)) {
    return(NULL)
  }
  # brmsframe prior objects: deparse to a stable character form.
  tryCatch(
    paste(capture.output(print(prior)), collapse = "\n"),
    error = function(e) NA_character_
  )
}

#' Normalize a brms stanvars object to a hashable representation
#'
#' @keywords internal
normalize_stanvars <- function(stanvars) {
  if (is.null(stanvars)) {
    return(NULL)
  }
  tryCatch(
    paste(capture.output(print(stanvars)), collapse = "\n"),
    error = function(e) NA_character_
  )
}

# ============================================================================
# Template data
# ============================================================================

#' Generate template data for prefit compilation
#'
#' Produces ONE small dataset using the configured data generator with a
#' throwaway RNG stream. The dataset is discarded after compilation; its only
#' purpose is to fix the Stan data structure (variable names and types) so the
#' compiled binary matches the real task data.
#'
#' Calls `data_generator(data_spec, task_ctx)` and returns its `$train` data
#' frame. `task_ctx$template` is TRUE so generators can recognize this
#' structural template call.
#'
#' @param data_generator The simulation data generator function.
#' @param data_spec A named list (one row of `data_grid` as a list).
#'
#' @return A data.frame suitable for `brms::brm(data =)`.
#'
#' @keywords internal
generate_template_data <- function(
  data_generator,
  data_spec
) {
  task_ctx <- list(
    task_id = "model_bank_template",
    data_idx = 1L,
    fit_idx = 1L,
    rep_idx = 0L,
    template = TRUE
  )

  data_bundle <- data_generator(data_spec, task_ctx)

  if (is.null(data_bundle) || is.null(data_bundle$train)) {
    stop(bayesim_internal_error(paste(
      "Template data generation did not return a data_bundle with $train.",
      "The model bank needs template data with the same structure as real",
      "task data to compile the Stan model."
    )))
  }

  if (!is.data.frame(data_bundle$train)) {
    stop(bayesim_internal_error(paste(
      "Template data_bundle$train must be a data.frame, got:",
      paste(class(data_bundle$train), collapse = ", ")
    )))
  }

  data_bundle$train
}

# ============================================================================
# Model bank construction
# ============================================================================

#' Resolve one model-grid row into canonical components
#'
#' @param fit_grid A model specification data frame.
#' @param i Scalar row index.
#' @return A named list with formula, family, prior, and stanvars.
#' @keywords internal
model_spec_from_grid_row <- function(fit_grid, i) {
  row <- as.list(fit_grid[i, , drop = FALSE])
  if (!"formula" %in% names(row) || is.null(row$formula)) {
    stop(bayesim_config_error(
      "fit_grid row " %+%
        i %+%
        " has no 'formula' column; every fit_grid row must specify a formula."
    ))
  }

  unwrap <- function(value, preserve_data_frame = FALSE) {
    if (
      length(value) == 1L &&
        is.list(value) &&
        !(preserve_data_frame && inherits(value, "data.frame"))
    ) {
      value[[1]]
    } else {
      value
    }
  }

  list(
    formula = unwrap(row$formula),
    family = unwrap(row$family),
    prior = unwrap(row$prior, preserve_data_frame = TRUE),
    stanvars = unwrap(row$stanvars)
  )
}

# Levels signature of categorical columns in a data frame.
#
# brms maps factor (and character) predictor values onto dummy codes via the
# column's LEVELS; under `update(recompile = FALSE)` the compiled binary
# reuses that mapping blindly. A task whose categorical predictors have
# different levels (or a different level ORDER, which silently reorders the
# coefficients) must therefore be treated as structurally incompatible with
# the prefit, even though `make_standata` field names and K are unchanged.
#
# Returns a named list of level vectors (factors: `levels()`; character:
# first-occurrence `unique()`, since brms factors them) or NULL when the data
# has no categorical columns. Name order and level order are significant.
.data_levels_sig <- function(data) {
  if (is.null(data) || !is.data.frame(data)) {
    return(NULL)
  }
  is_cat <- vapply(
    data,
    function(col) is.factor(col) || is.character(col),
    logical(1)
  )
  if (!any(is_cat)) {
    return(NULL)
  }
  lapply(data[is_cat], function(col) {
    if (is.factor(col)) {
      levels(col)
    } else {
      unique(col)
    }
  })
}

#' Build the model bank for a BrmsFitter
#'
#' For each DISTINCT row of `fit_grid` (deduped by [model_spec_hash()]), compiles
#' a prefit via `brms::brm(chains = 0)` against generator-supplied template
#' data. Returns a named list of prefit objects keyed by spec hash.
#'
#' When `fitter@precompile` is FALSE, returns NULL immediately (the fitter falls
#' back to a fresh `brms::brm()` per task).
#'
#' A compile failure is fatal: it raises a [bayesim_internal_error()] since the
#' simulation cannot proceed without a compilable model.
#'
#' @param fitter A BrmsFitter S7 object.
#' @param fit_grid A data.frame of model fitting specifications. Each row's
#'   `formula`, `family`, `prior`, `stanvars` columns (list-columns) define a
#'   model spec.
#' @param data_generator The simulation data generator function.
#' @param data_spec_template A named list (one row of `data_grid` as a list)
#'   used to generate template data.
#' @param result_path NULL or character path. When non-NULL,
#'   `options(cmdstanr_write_stan_file_dir)` is set to
#'   `file.path(result_path, "stan_binaries")` so compiled binaries persist and
#'   are shared across controller and local daemons.
#'
#' @return A named list of `brmsfit` prefit objects keyed by
#'   [model_spec_hash()], or NULL when `precompile` is FALSE.
#'
#' @keywords internal
build_model_bank <- function(
  fitter,
  fit_grid,
  data_generator,
  data_spec_template,
  result_path = NULL
) {
  if (!isTRUE(fitter@precompile)) {
    return(NULL)
  }

  backend <- fitter@backend

  # Generate template data once (same data structure for all specs; the spec
  # determines Stan code, not the data values).
  template_data <- generate_template_data(
    data_generator = data_generator,
    data_spec = data_spec_template
  )

  # Configure persistent cmdstanr binary output dir so binaries survive across
  # processes (controller + daemons share the filesystem).
  if (!is.null(result_path)) {
    bin_dir <- file.path(result_path, "stan_binaries")
    dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
    withr::local_options(list(cmdstanr_write_stan_file_dir = bin_dir))
  }

  # Iterate distinct specs. Dedupe by hash so identical specs compile once.
  model_bank <- list()

  for (i in seq_len(nrow(fit_grid))) {
    spec <- model_spec_from_grid_row(fit_grid, i)
    formula <- spec$formula
    family <- spec$family
    prior <- spec$prior
    stanvars <- spec$stanvars

    if (is.null(prior) || length(prior) == 0L) {
      .warn_once(
        "model_bank_default_prior",
        c(
          "Precompiled brms models should use explicit priors.",
          i = "Some brms default priors depend on the template data and remain embedded when the compiled model is reused.",
          i = "Set an explicit {.code prior} in every model-grid row, or use {.code BrmsFitter(precompile = FALSE)}."
        )
      )
    }

    spec_hash <- model_spec_hash(
      formula = formula,
      family = family,
      prior = prior,
      stanvars = stanvars,
      backend = backend
    )

    if (spec_hash %in% names(model_bank)) {
      next # already compiled an identical spec
    }

    formula_resolved <- resolve_formula(formula)
    family_resolved <- resolve_family(family)

    cli::cli_alert_info(
      "Compiling model spec {length(model_bank) + 1L} (hash {substr(spec_hash, 1, 8)}): {.code {deparse(formula_resolved$formula)[1]}}"
    )

    prefit <- tryCatch(
      brms::brm(
        formula = formula_resolved,
        data = template_data,
        family = family_resolved,
        prior = prior,
        stanvars = stanvars,
        backend = backend,
        chains = 0L,
        refresh = 0L,
        silent = 2L,
        threads = stan_threads_arg(fitter@stan_args)
      ),
      error = function(e) e
    )

    if (inherits(prefit, "error") || inherits(prefit, "simpleError")) {
      stop(bayesim_internal_error(paste(
        "Model bank compilation failed for spec",
        sprintf("(hash %s):", substr(spec_hash, 1, 8)),
        conditionMessage(prefit)
      )))
    }

    # F6: precompute the prefit-side Stan data STRUCTURE signature once here,
    # so update_prefit() does not recompute brms::make_standata() per task
    # (200x for a 10k-task run). Same struct_sig definition as in update_prefit.
    prefit_struct <- tryCatch(
      brms::make_standata(
        formula = formula_resolved,
        data = prefit$data,
        family = family_resolved
      ),
      error = function(e) NULL
    )
    struct_sig <- if (!is.null(prefit_struct)) {
      list(
        fields = names(prefit_struct),
        K = prefit_struct$K,
        X_ncol = ncol(prefit_struct$X),
        levels = .data_levels_sig(prefit$data)
      )
    } else {
      NULL
    }

    model_bank[[spec_hash]] <- list(prefit = prefit, struct_sig = struct_sig)
  }

  cli::cli_alert_info(
    "Model bank built: {length(model_bank)} distinct model spec(s) compiled"
  )

  model_bank
}

#' Extract the `threads` argument from stan_args, if present
#'
#' @keywords internal
stan_threads_arg <- function(stan_args) {
  if (is.null(stan_args)) {
    return(NULL)
  }
  if (!is.null(stan_args$threads)) {
    return(stan_args$threads)
  }
  NULL
}

# ============================================================================
# Prefit lookup
# ============================================================================

#' Look up a prefit in the model bank
#'
#' Computes the spec hash and returns the matching prefit object, or NULL if
#' the spec is not in the bank (caller falls back to a fresh compile).
#'
#' @param model_bank A named list of prefit objects (from [build_model_bank()]),
#'   or NULL.
#' @param formula A formula, brmsformula, or string.
#' @param family A brms family object, string, or NULL.
#' @param prior A brms prior object or NULL.
#' @param stanvars A brms stanvars object or NULL.
#' @param backend Character scalar.
#'
#' @return A `brmsfit` prefit object, or NULL.
#'
#' @keywords internal
lookup_prefit <- function(
  model_bank,
  formula,
  family,
  prior,
  stanvars,
  backend
) {
  if (is.null(model_bank) || length(model_bank) == 0L) {
    return(NULL)
  }

  spec_hash <- model_spec_hash(
    formula = formula,
    family = family,
    prior = prior,
    stanvars = stanvars,
    backend = backend
  )

  model_bank[[spec_hash]] %||% NULL
}
