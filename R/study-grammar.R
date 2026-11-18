# Study grammar --------------------------------------------------------------
#
# A study declares conditions, a data generator, methods, measurements,
# comparisons and what to keep. The declaration is a plain list: nothing runs
# until run_study(). plan_study() checks the declaration and captures its
# functions; the engine in study-run.R executes it.

#' Declare a simulation study
#'
#' `study()` starts a study declaration. Add methods with [with_methods()],
#' per-fit measurements with [with_measure()], comparisons with
#' [with_comparison()] and artifacts to keep with [with_retention()]. Check
#' the declaration with [plan_study()] and run it with [run_study()].
#'
#' Each dataset is identified by the seed, the condition row, the replicate
#' number, the value of `prepare` and the generator. All methods see the same
#' dataset for a condition and replicate.
#'
#' @param name Study name, for printing and reports. The functions of the
#'   study do not see it, so studies that differ only in name share results
#'   in one store.
#' @param conditions A data frame with one row per condition. It must have a
#'   unique, nonempty character column `condition_id`. Atomic columns appear
#'   in results with the prefix `condition_`. A column named `id` is not
#'   allowed, because it would appear as `condition_id`.
#' @param generate `function(condition, context)` that returns one dataset.
#'   `condition` is the condition row as a list. `context` holds
#'   `condition`, `dataset_id`, `replicate` and `prepared`.
#' @param prepare Optional `function(condition, context)` that runs once per
#'   condition. Its value is available as `context$prepared` to every other
#'   function.
#' @param version Version of `generate` and `prepare`. Change it when they
#'   depend on something that bayesim cannot see, such as a file.
#'
#' @return A `bayesim_study_declaration`.
#' @seealso [study_method()], [run_study()]
#' @export
#' @examples
#' s <- study(
#'   "normal mean",
#'   data.frame(condition_id = c("small", "large"), n = c(10, 100)),
#'   prior_draws(
#'     prior = function(condition) c(mu = rnorm(1)),
#'     simulate = function(truth, condition) {
#'       X <- matrix(1, condition$n, 1, dimnames = list(NULL, "mu"))
#'       list(X = X, y = rnorm(condition$n, truth[["mu"]]))
#'     }
#'   )
#' )
#' s <- with_methods(s, conjugate = conjugate_lm())
#' s <- with_measure(s, "coverage", measure_coverage(0.9))
#' run <- run_study(s, replicates = 20)
#' summarise_study(run)
study <- function(name, conditions, generate, prepare = NULL, version = "1") {
  .sg_name(name, "name")
  if (
    !is.data.frame(conditions) ||
      !nrow(conditions) ||
      !is.character(conditions$condition_id) ||
      anyNA(conditions$condition_id) ||
      !all(nzchar(conditions$condition_id)) ||
      anyDuplicated(conditions$condition_id)
  ) {
    cli::cli_abort(
      "{.arg conditions} needs a unique, nonempty character {.field condition_id} column."
    )
  }
  .sg_function(generate, "generate")
  if (!is.null(prepare)) {
    .sg_function(prepare, "prepare")
  }
  .sg_name(version, "version")
  structure(
    list(
      name = name,
      conditions = conditions,
      generate = generate,
      prepare = prepare,
      version = version,
      methods = list(),
      measures = list(),
      comparisons = list(),
      retention = character(),
      compress = "gzip"
    ),
    class = "bayesim_study_declaration"
  )
}

#' Declare a study method
#'
#' A method fits a model to one dataset and names the artifacts that
#' measurements and comparisons may use.
#'
#' The fit runs on its own random stream. The stream depends on the seed, the
#' dataset and the method ID, so two methods never share random numbers and
#' a method gets the same numbers in every run. Functions can run in other R
#' processes: bayesim copies the values they use from outside packages, and
#' [plan_study()] stops if a value cannot be copied or is too large.
#'
#' @param fit `function(data, context)` that returns a fitted model.
#'   `context` adds `method_id` and `settings` to the generator's context. An
#'   error in `fit` is recorded as a failed attempt; the study continues.
#' @param extract A named list of `function(fit, data, context)`. Each one
#'   computes an artifact from the fit. The names `fit` and `data` are
#'   reserved for the fit and the dataset.
#' @param settings A named list. It is available as `context$settings`.
#'   Atomic settings of length one appear in results with the prefix
#'   `method_`. A setting named `id` is not allowed, because it would appear
#'   as `method_id`.
#' @param version Version of `fit` and `extract`. Change it when they depend
#'   on something that bayesim cannot see, such as a file or a compiled
#'   model.
#'
#' @return A `bayesim_study_method`.
#' @export
study_method <- function(
  fit,
  extract = list(),
  settings = list(),
  version = "1"
) {
  .sg_function(fit, "fit")
  .sg_named(extract, "extract")
  for (nm in names(extract)) {
    .sg_function(extract[[nm]], paste0("extract$", nm))
  }
  if (any(names(extract) %in% c("fit", "data"))) {
    cli::cli_abort(
      "{.code fit} and {.code data} are built-in artifacts. Use other extractor names."
    )
  }
  if (!is.list(settings)) {
    cli::cli_abort("{.arg settings} must be a list.")
  }
  .sg_named(settings, "settings")
  .sg_name(version, "version")
  structure(
    list(fit = fit, extract = extract, settings = settings, version = version),
    class = "bayesim_study_method"
  )
}

#' Add methods to a study
#'
#' @param study A [study()] declaration.
#' @param ... Named [study_method()] objects. The names are the method IDs.
#' @return The declaration with the methods added.
#' @export
with_methods <- function(study, ...) {
  .sg_study(study)
  methods <- list(...)
  .sg_named(methods, "methods")
  if (
    !length(methods) ||
      !all(vapply(methods, inherits, logical(1), "bayesim_study_method"))
  ) {
    cli::cli_abort("Supply named {.fn study_method} objects.")
  }
  if (any(names(methods) %in% names(study$methods))) {
    cli::cli_abort("Method IDs already exist.")
  }
  study$methods <- c(study$methods, methods)
  study
}

#' Add a measurement to a study
#'
#' A measurement runs once for each successful fit. It returns a data frame,
#' usually with columns `target` and `value`. bayesim adds the condition,
#' dataset, replicate and method columns and the column `measurement`.
#'
#' @param study A [study()] declaration.
#' @param name Measurement name.
#' @param compute `function(artifacts, context)`. `artifacts` is a named
#'   list with the artifacts in `needs`. Each measurement has its own random
#'   stream.
#' @param needs Names of the artifacts that `compute` uses: `fit`, `data` or
#'   an extractor name. Built-in measurements such as [measure_coverage()] supply
#'   their own default.
#' @param version Version of `compute`.
#' @return The declaration with the measurement added.
#' @export
with_measure <- function(
  study,
  name,
  compute,
  needs = attr(compute, "needs") %||% character(),
  version = "1"
) {
  .sg_study(study)
  .sg_name(name, "name")
  .sg_function(compute, "compute")
  .sg_needs(needs, "needs")
  .sg_name(version, "version")
  if (name %in% names(study$measures)) {
    cli::cli_abort("Measurement name already exists.")
  }
  study$measures[[name]] <- list(
    compute = compute,
    needs = unique(needs),
    version = version,
    summary = attr(compute, "summary") %||% "mean"
  )
  study
}

#' Add a comparison to a study
#'
#' A comparison runs once for each group of attempts defined by `by`.
#'
#' A comparison without `needs` uses only measurements. It runs after the
#' fits, and [assess_study()] runs it again on the measurements that an
#' analysis policy keeps. A comparison with `needs` uses artifacts. It runs
#' as soon as every method has been fitted to a dataset, so `by` must
#' include `dataset_id`. Its result is stored, and artifacts that only it
#' needs are deleted afterwards.
#'
#' @param study A [study()] declaration.
#' @param name Comparison name.
#' @param compute `function(results, artifacts, context)`. `results` holds
#'   the measurement rows of the group. `artifacts` is a list by method ID of
#'   the artifacts in `needs`, for successful fits only. `context$attempts`
#'   lists every attempt in the group, including failed ones, and
#'   `context$group` holds the values of `by`.
#' @param by Result columns that define the groups.
#' @param needs Names of the artifacts that `compute` uses.
#' @param version Version of `compute`.
#' @return The declaration with the comparison added.
#' @export
with_comparison <- function(
  study,
  name,
  compute,
  by = "dataset_id",
  needs = character(),
  version = "1"
) {
  .sg_study(study)
  .sg_name(name, "name")
  .sg_function(compute, "compute")
  .sg_needs(needs, "needs")
  .sg_needs(by, "by")
  .sg_name(version, "version")
  if (name %in% names(study$comparisons)) {
    cli::cli_abort("Comparison name already exists.")
  }
  if (length(needs) && !"dataset_id" %in% by) {
    cli::cli_abort(
      "Artifact-consuming comparisons must group within dataset_id."
    )
  }
  study$comparisons[[name]] <- list(
    compute = compute,
    by = by,
    needs = unique(needs),
    version = version
  )
  study
}

#' Keep artifacts after a fit
#'
#' Kept artifacts let later runs add measurements without fitting again.
#' Keep `fit` to store the native fit object. Keep `data` to store each
#' dataset once instead of generating it again. Artifacts that are not kept
#' are deleted when no measurement or comparison needs them any more.
#'
#' @param study A [study()] declaration.
#' @param artifacts Names of artifacts to keep.
#' @param compress Compression for stored artifacts: `"gzip"`, `"bzip2"`
#'   or `"xz"`.
#' @return The declaration with the retention set.
#' @export
with_retention <- function(study, artifacts = character(), compress = "gzip") {
  .sg_study(study)
  .sg_needs(artifacts, "artifacts")
  if (
    !is.character(compress) ||
      length(compress) != 1L ||
      !compress %in% c("gzip", "bzip2", "xz")
  ) {
    cli::cli_abort("Use gzip, bzip2 or xz compression.")
  }
  study$retention <- unique(artifacts)
  study$compress <- compress
  study
}

#' Check a study and describe the work
#'
#' `plan_study()` runs no data generation, fits or measurements. It checks
#' that every method can supply the artifacts that measurements, comparisons
#' and retention need, that comparison groups use columns that results will
#' have, and that every function can be sent to another R process.
#'
#' @param study A [study()] declaration.
#' @param replicates Number of replicates per condition.
#' @param seed Study seed, a whole number.
#' @param path Optional store directory. When it is given, the plan also
#'   counts the fits that are done, failed and still to run.
#'
#' @return A `bayesim_study_plan`.
#' @export
plan_study <- function(study, replicates, seed = 1L, path = NULL) {
  if (is.null(path)) {
    plan <- .sg_plan(study, replicates, seed)
    attr(plan, "capture") <- NULL
  } else {
    setup <- .sg_setup(study, replicates, seed)
    plan <- setup$plan
    st <- .store_open(path)
    jobs <- .sg_fit_states(setup, .store_index(st, "fits")$fits, retry = FALSE)
    plan$store <- c(
      done = sum(jobs$state == "done"),
      failed = sum(jobs$state == "failed"),
      measure = sum(jobs$state == "measure"),
      fit = sum(jobs$state == "fit")
    )
  }
  plan
}

#' @export
print.bayesim_study_plan <- function(x, ...) {
  cat(
    x$study,
    "\n",
    x$conditions,
    " conditions x ",
    x$replicates,
    " replicates = ",
    format(x$datasets, scientific = FALSE),
    " datasets\n",
    length(x$methods),
    " methods; ",
    format(x$fits, scientific = FALSE),
    " fits\n",
    sep = ""
  )
  if (x$prepare) {
    cat("Prepare once per condition.\n")
  }
  print(x$dependencies, row.names = FALSE)
  cat(
    "Retain: ",
    if (length(x$retained)) {
      paste(x$retained, collapse = ", ")
    } else {
      "measurements only"
    },
    "\n",
    sep = ""
  )
  cat(
    "Discard after use: ",
    paste(x$discarded, collapse = ", "),
    "\n",
    sep = ""
  )
  if (!is.null(x$store)) {
    cat(
      "Store: ",
      x$store[["done"]],
      " done, ",
      x$store[["failed"]],
      " failed, ",
      x$store[["measure"]],
      " to measure again, ",
      x$store[["fit"]],
      " to fit.\n",
      sep = ""
    )
  } else {
    cat("Plan only: no data, models or store records have been evaluated.\n")
  }
  invisible(x)
}

#' @export
print.bayesim_study_declaration <- function(x, ...) {
  cat(
    "<study> ",
    x$name,
    ": ",
    nrow(x$conditions),
    " conditions, ",
    length(x$methods),
    " methods, ",
    length(x$measures),
    " measurements, ",
    length(x$comparisons),
    " comparisons\n",
    sep = ""
  )
  invisible(x)
}

# Validation -----------------------------------------------------------------

.sg_plan <- function(study, replicates, seed) {
  .sg_study(study)
  .sg_count(replicates, "replicates")
  .sg_count(seed, "seed", zero = TRUE)
  if (!length(study$methods)) {
    cli::cli_abort("Add at least one method.")
  }
  .sg_reserved_fields(study)
  needs <- .sg_required(study)
  by <- unique(unlist(lapply(study$comparisons, `[[`, "by"), use.names = FALSE))
  columns <- c(
    "condition_id",
    "dataset_id",
    "replicate",
    "method_id",
    names(.sg_condition_meta(study$conditions)),
    names(.sg_method_meta(study$methods))
  )
  unknown <- setdiff(by, columns)
  if (length(unknown)) {
    cli::cli_abort(
      "Unknown comparison grouping columns: {paste(unknown, collapse = ', ')}."
    )
  }
  for (id in names(study$methods)) {
    available <- c("data", "fit", names(study$methods[[id]]$extract))
    missing <- setdiff(needs, available)
    if (length(missing)) {
      cli::cli_abort(
        "Method {.val {id}} cannot supply: {paste(missing, collapse = ', ')}."
      )
    }
  }
  # Capture now, so that functions which cannot travel fail before any work.
  captured <- .sg_capture(study)
  dependencies <- .sg_bind(c(
    lapply(names(study$measures), function(id) {
      data.frame(
        stage = "per fit",
        name = id,
        needs = paste(study$measures[[id]]$needs, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }),
    lapply(names(study$comparisons), function(id) {
      x <- study$comparisons[[id]]
      data.frame(
        stage = if (length(x$needs)) {
          "within dataset, before eviction"
        } else {
          "saved measurements"
        },
        name = id,
        needs = paste(x$needs, collapse = ", "),
        stringsAsFactors = FALSE
      )
    })
  ), order = FALSE)
  structure(
    list(
      study = study$name,
      conditions = nrow(study$conditions),
      replicates = as.integer(replicates),
      datasets = nrow(study$conditions) * as.double(replicates),
      methods = names(study$methods),
      fits = nrow(study$conditions) *
        as.double(replicates) *
        length(study$methods),
      prepare = !is.null(study$prepare),
      required = needs,
      retained = study$retention,
      discarded = setdiff(needs, study$retention),
      dependencies = dependencies,
      seed = as.integer(seed)
    ),
    class = "bayesim_study_plan",
    capture = captured
  )
}

# Issue #110: a field named `id` would appear as `condition_id` or
# `method_id` and overwrite the identity column. Refuse it instead.
.sg_reserved_fields <- function(study) {
  if ("id" %in% names(study$conditions)) {
    cli::cli_abort(c(
      "Condition field {.field id} would appear in results as {.field condition_id}, which identifies the condition.",
      i = "Rename the field."
    ))
  }
  for (method in names(study$methods)) {
    if ("id" %in% names(study$methods[[method]]$settings)) {
      cli::cli_abort(c(
        "Method {.val {method}} has a setting {.field id}, which would appear in results as {.field method_id}.",
        i = "Rename the setting."
      ))
    }
  }
}

.sg_name <- function(x, arg = "name") {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    cli::cli_abort("{.arg {arg}} must be one nonempty name or version.")
  }
}
.sg_function <- function(x, arg = "x") {
  if (!is.function(x)) {
    cli::cli_abort("{.arg {arg}} must be a function.")
  }
}
.sg_needs <- function(x, arg = "needs") {
  if (!is.character(x) || anyNA(x) || !all(nzchar(x)) || anyDuplicated(x)) {
    cli::cli_abort("{.arg {arg}} must hold unique, nonempty names.")
  }
}
.sg_named <- function(x, arg = "x") {
  if (
    length(x) &&
      (is.null(names(x)) ||
        anyNA(names(x)) ||
        !all(nzchar(names(x))) ||
        anyDuplicated(names(x)))
  ) {
    cli::cli_abort(
      "Every entry of {.arg {arg}} must have a unique, nonempty name."
    )
  }
}
.sg_study <- function(x) {
  if (!inherits(x, "bayesim_study_declaration")) {
    cli::cli_abort("Expected a {.fn study} declaration.")
  }
}
.sg_count <- function(x, name, zero = FALSE) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x != floor(x) ||
      x < (if (zero) 0 else 1) ||
      x > .Machine$integer.max
  ) {
    cli::cli_abort(
      "{.arg {name}} must be a whole number in the supported integer range."
    )
  }
}

# Artifacts that a fit must produce.
.sg_required <- function(study) {
  unique(c(
    study$retention,
    unlist(
      lapply(c(study$measures, study$comparisons), `[[`, "needs"),
      use.names = FALSE
    )
  ))
}

# Result metadata ------------------------------------------------------------
#
# Atomic, non-matrix condition columns become `condition_<name>`; atomic
# settings of length one become `method_<name>`. Identity columns come first
# and nothing can overwrite them (#110).

.sg_condition_meta <- function(conditions) {
  keep <- vapply(
    conditions,
    function(x) is.atomic(x) && is.null(dim(x)),
    logical(1)
  )
  keep[names(conditions) %in% c("condition_id", "id")] <- FALSE
  out <- lapply(conditions[keep], function(x) {
    if (inherits(x, "AsIs")) unclass(x) else x
  })
  names(out) <- .sg_prefix("condition_", names(out))
  .sg_frame(out, nrow(conditions))
}

.sg_method_meta <- function(methods) {
  values <- lapply(methods, function(m) {
    Filter(
      function(x) is.atomic(x) && length(x) == 1L,
      m$settings[setdiff(names(m$settings), "id")]
    )
  })
  columns <- unique(unlist(lapply(values, names), use.names = FALSE))
  out <- lapply(columns, function(nm) {
    do.call(vctrs::vec_c, unname(lapply(values, function(v) v[[nm]] %||% NA)))
  })
  names(out) <- .sg_prefix("method_", columns)
  .sg_frame(out[order(names(out))], length(methods))
}

# paste0() returns the prefix alone for an empty vector; this returns nothing.
.sg_prefix <- function(prefix, x) {
  if (length(x)) paste0(prefix, x) else character()
}

# A data frame from a list of columns, also with no columns.
.sg_frame <- function(columns, n) {
  columns <- as.list(columns)
  attr(columns, "row.names") <- .set_row_names(n)
  class(columns) <- "data.frame"
  columns
}

# Bind data frames that may have different columns. Missing columns become
# NA, or NULL in list columns. Identity columns come first and the others are
# sorted, so the column order does not depend on which rows exist.
.sg_identities <- c(
  "condition_id",
  "dataset_id",
  "replicate",
  "method_id",
  "measurement",
  "comparison"
)

.sg_bind <- function(rows, order = TRUE, keep_empty = FALSE) {
  rows <- Filter(is.data.frame, rows)
  if (!keep_empty) {
    rows <- Filter(function(x) nrow(x) > 0L, rows)
  }
  if (!length(rows)) {
    return(data.frame())
  }
  out <- if (length(rows) == 1L) {
    rows[[1L]]
  } else {
    do.call(vctrs::vec_rbind, unname(rows))
  }
  out <- .sg_frame(as.list(out), vctrs::vec_size(out))
  if (order) {
    out <- .sg_order_columns(out)
  }
  out
}

.sg_order_columns <- function(x) {
  columns <- names(x)
  x[c(
    intersect(.sg_identities, columns),
    sort(setdiff(columns, .sg_identities))
  )]
}
