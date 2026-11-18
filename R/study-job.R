# Study jobs -------------------------------------------------------------------
#
# A job is one unit of work that a backend runs: prepare a condition, fit one
# method to one dataset (or measure an existing fit again), or run the
# artifact comparisons of one dataset. Workers read inputs from the store and
# write their records to it. They return only a short status, so the
# controller never receives fits, artifacts or measurement rows.
#
# The worker is a closure over the run specification. A backend sends it to
# each process once; each job is then a short list of indices and keys.

.make_worker <- function(spec) {
  force(spec)
  cache <- new.env(parent = emptyenv())
  function(job) .run_job(job, spec, cache)
}

# Run one job. An error in `fit` becomes a failed record. Any other error is
# returned with the stage, so that the controller can stop the run and say
# where the error happened. Warnings are recorded with the attempt.
.run_job <- function(job, spec, cache) {
  # With a state in place, each stream swaps .Random.seed and nothing else.
  # The caller's state is restored by run_study().
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    stats::runif(1L)
  }
  cache$stage <- job$kind
  cache$warnings <- character()
  started <- proc.time()[["elapsed"]]
  out <- tryCatch(
    withCallingHandlers(
      switch(
        job$kind,
        prepare = .job_prepare(job, spec, cache),
        fit = .job_fit(job, spec, cache),
        compare = .job_compare(job, spec, cache)
      ),
      warning = function(w) {
        cache$warnings <- c(cache$warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(
        error = list(
          message = conditionMessage(e),
          stage = cache$stage,
          condition_id = spec$condition_ids[[job$cidx]],
          replicate = job$replicate,
          method_id = if (job$kind == "fit") spec$method_ids[[job$midx]]
        )
      )
    }
  )
  out$kind <- job$kind
  out$seconds <- proc.time()[["elapsed"]] - started
  out$warnings <- length(cache$warnings)
  out
}

# Inputs ----------------------------------------------------------------------

.job_store <- function(spec) .store_handle(spec$root)

.job_prepared <- function(spec, cache, cidx) {
  if (is.null(spec$cap$prepare)) {
    return(NULL)
  }
  key <- spec$prepare_keys[[cidx]]
  cache$prepared <- cache$prepared %||% list()
  if (!is.null(cache$prepared[[key]])) {
    return(cache$prepared[[key]]$value)
  }
  file <- .store_file(.job_store(spec), "prepared", key)
  if (!.store_has(file)) {
    cli::cli_abort(
      "The prepared value for condition {.val {spec$condition_ids[[cidx]]}} is missing from the store."
    )
  }
  value <- .read_payload(file)
  # Keep the last few conditions; a worker usually sees one at a time.
  cache$prepared <- utils::tail(c(cache$prepared, stats::setNames(list(value), key)), 8L)
  value$value
}

.job_context <- function(job, spec, cache) {
  list(
    study = spec$cap$name,
    condition = spec$conditions[[job$cidx]],
    dataset_id = job$dataset_id,
    replicate = job$replicate,
    prepared = .job_prepared(spec, cache, job$cidx)
  )
}

.job_method_context <- function(context, spec, midx) {
  context$method_id <- spec$method_ids[[midx]]
  context$settings <- spec$cap$methods[[midx]]$settings
  context
}

# The dataset: from this process's cache, from the store when it is retained,
# or generated again on its own stream.
.job_data <- function(job, context, spec, cache) {
  if (identical(cache$data_id, job$dataset_id)) {
    return(cache$data)
  }
  st <- .job_store(spec)
  file <- .store_file(st, "datasets", job$dataset_id)
  data <- if (.store_has(file)) {
    .read_payload(file)$data
  } else {
    cache$stage <- "generate"
    value <- .invoke(
      c("data", job$dataset_id),
      spec$cap$generate,
      condition = context$condition,
      context = context
    )
    if ("data" %in% spec$retention) {
      .write_payload(list(data = value), file, compress = spec$compress)
    }
    value
  }
  cache$data_id <- job$dataset_id
  cache$data <- data
  data
}

# Prepare ----------------------------------------------------------------------

.job_prepare <- function(job, spec, cache) {
  condition <- spec$conditions[[job$cidx]]
  key <- spec$prepare_keys[[job$cidx]]
  value <- .invoke(
    c("prepare", key),
    spec$cap$prepare,
    condition = condition,
    context = list(study = spec$cap$name, condition = condition)
  )
  .write_payload(
    list(value = value),
    .store_file(.job_store(spec), "prepared", key),
    compress = spec$compress
  )
  list(status = "success")
}

# Fit and measure ------------------------------------------------------------------

.job_fit <- function(job, spec, cache) {
  st <- .job_store(spec)
  method <- spec$cap$methods[[job$midx]]
  method_id <- spec$method_ids[[job$midx]]
  base <- .job_context(job, spec, cache)
  context <- .job_method_context(base, spec, job$midx)
  data <- NULL
  get_data <- function() {
    if (is.null(data)) {
      data <<- list(.job_data(job, base, spec, cache))
    }
    data[[1L]]
  }

  if (identical(job$mode, "remeasure")) {
    old <- .read_payload(.store_file(st, "fits", job$old))
    if (!identical(old$extract_id, method$extract_id)) {
      cli::cli_abort(c(
        "Extractor recipe changed for cached method {.val {method_id}}.",
        i = "Change the method version to fit again."
      ))
    }
    attempt <- old$attempt
    cached <- old$measurements
    artifacts <- .read_if(.store_file(st, "artifacts", job$fit_key)) %||%
      list()
    missing <- setdiff(spec$retention, c("data", names(artifacts)))
    if (length(missing)) {
      cli::cli_abort(c(
        "Requested retention needs discarded artifacts: {paste(missing, collapse = ', ')}.",
        i = "Use a new store path to fit again with these artifacts kept."
      ))
    }
    kept <- NULL
  } else {
    get_data()
    cache$stage <- "fit"
    fit_started <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      list(
        value = .invoke(
          c("fit", job$dataset_id, method_id),
          method$fit,
          data = data[[1L]],
          context = context
        )
      ),
      error = function(e) list(error = conditionMessage(e))
    )
    attempt <- list(
      status = if (is.null(fit$error)) "success" else "failed",
      error = fit$error %||% NA_character_,
      elapsed = proc.time()[["elapsed"]] - fit_started,
      warnings = .job_warnings(cache)
    )
    failed_file <- .store_file(st, "fits", job$fit_key, ".failed.rds")
    if (!is.null(fit$error)) {
      .write_payload(
        list(attempt = attempt, extract_id = method$extract_id),
        failed_file
      )
      return(list(status = "failed", record = paste0(job$fit_key, ".failed")))
    }
    artifacts <- list()
    cache$stage <- "extract"
    for (name in spec$required[spec$required != "data"]) {
      artifacts[name] <- list(
        if (name == "fit") {
          fit$value
        } else {
          .invoke(
            c("extract", job$dataset_id, method_id, name),
            method$extract[[name]],
            fit = fit$value,
            data = data[[1L]],
            context = context
          )
        }
      )
    }
    rm(fit)
    cached <- list()
    kept <- unique(c(spec$retention, spec$comparison_needs))
    kept <- kept[kept %in% names(artifacts)]
  }

  measurements <- list()
  for (name in names(spec$cap$measures)) {
    measure <- spec$cap$measures[[name]]
    old <- cached[[name]]
    if (!is.null(old) && identical(old$id, measure$id)) {
      measurements[[name]] <- old
      next
    }
    cache$stage <- "measure"
    if ("data" %in% measure$needs) {
      artifacts["data"] <- list(get_data())
    }
    missing <- measure$needs[!measure$needs %in% names(artifacts)]
    if (length(missing)) {
      cli::cli_abort(c(
        "Measurement {.val {name}} needs discarded artifacts: {paste(missing, collapse = ', ')}.",
        i = "Use a new store path to fit again with these artifacts kept."
      ))
    }
    rows <- .invoke(
      c("measure", job$dataset_id, method_id, name),
      measure$compute,
      artifacts = artifacts[measure$needs],
      context = context
    )
    measurements[[name]] <- list(
      id = measure$id,
      rows = .job_check_rows(rows, name, spec)
    )
  }
  if (!is.null(kept)) {
    if (length(kept)) {
      .write_payload(
        artifacts[kept],
        .store_file(st, "artifacts", job$fit_key),
        compress = spec$compress
      )
    }
  }
  attempt$warnings <- .job_warnings(cache, attempt$warnings)
  record <- paste0(job$fit_key, "-", job$mhash)
  .write_payload(
    list(
      attempt = attempt,
      measurements = measurements,
      extract_id = method$extract_id
    ),
    .store_file(st, "fits", record)
  )
  unlink(.store_file(st, "fits", job$fit_key, ".failed.rds"))
  list(status = "success", record = record)
}

.job_warnings <- function(cache, old = NA_character_) {
  if (is.na(old) && !length(cache$warnings)) {
    return(NA_character_)
  }
  w <- unique(c(if (!is.na(old)) old, cache$warnings))
  if (length(w)) paste(w, collapse = "\n") else NA_character_
}

.job_check_rows <- function(rows, name, spec) {
  if (!is.data.frame(rows)) {
    cli::cli_abort("Measurement {.val {name}} must return a data frame.")
  }
  if (any(names(rows) %in% c(spec$meta_names, "measurement"))) {
    cli::cli_abort(
      "Measurement output uses reserved identity columns: {paste(intersect(names(rows), c(spec$meta_names, 'measurement')), collapse = ', ')}."
    )
  }
  rows <- .sg_frame(as.list(rows), nrow(rows))
  rows
}

# Artifact comparisons ------------------------------------------------------------

# Run the artifact comparisons of one dataset after all its methods are
# fitted, store each result, then delete artifacts that only comparisons
# needed.
.job_compare <- function(job, spec, cache) {
  st <- .job_store(spec)
  records <- lapply(job$records, function(r) {
    .read_payload(.store_file(st, "fits", r))
  })
  n <- length(job$records)
  status <- vapply(records, function(r) r$attempt$status, character(1))
  ix <- list(
    didx = rep(1L, n),
    dataset_id = rep(job$dataset_id, n),
    cidx = rep(job$cidx, n),
    replicate = rep(job$replicate, n),
    midx = job$midx,
    fit_key = job$fit_keys
  )
  attempts <- .sg_attempt_frame(spec, ix, lapply(records, `[[`, "attempt"))
  results <- .sg_measure_frame(
    spec,
    ix,
    lapply(records, function(r) r$measurements %||% list())
  )
  context <- .job_context(job, spec, cache)
  bundles <- list()
  for (i in which(status == "success")) {
    a <- .read_if(.store_file(st, "artifacts", job$fit_keys[[i]])) %||% list()
    if ("data" %in% spec$comparison_needs) {
      a["data"] <- list(.job_data(job, context, spec, cache))
    }
    bundles[[spec$method_ids[[job$midx[[i]]]]]] <- a
  }
  cache$stage <- "compare"
  for (name in names(job$keys)) {
    file <- .store_file(st, "comparisons", job$keys[[name]])
    if (.store_has(file)) {
      next
    }
    value <- .sg_compare(
      spec$cap$comparisons[name],
      spec$cap$name,
      results,
      attempts,
      bundles = bundles,
      stream = c("compare", job$dataset_id)
    )
    .write_payload(list(value = value[[name]]), file)
  }
  # Evict artifacts that no later step needs.
  for (i in which(status == "success")) {
    file <- .store_file(st, "artifacts", job$fit_keys[[i]])
    if (!.store_has(file)) {
      next
    }
    a <- .read_payload(file)
    keep <- intersect(names(a), spec$retention)
    if (!length(keep)) {
      unlink(file)
    } else if (length(keep) < length(a)) {
      .write_payload(a[keep], file, compress = spec$compress)
    }
  }
  list(status = "success")
}

# Result tables ---------------------------------------------------------------------
#
# Records hold only what a fit produced. Identity and metadata columns are
# added here, from index vectors, so that the controller and comparison jobs
# build the same tables.

.sg_meta_frame <- function(spec, ix) {
  cbind(
    .sg_frame(
      list(
        condition_id = spec$condition_ids[ix$cidx],
        dataset_id = ix$dataset_id,
        replicate = ix$replicate,
        method_id = spec$method_ids[ix$midx]
      ),
      length(ix$cidx)
    ),
    spec$cond_meta[ix$cidx, , drop = FALSE],
    spec$method_meta[ix$midx, , drop = FALSE]
  )
}

.sg_attempt_frame <- function(spec, ix, attempts) {
  n <- length(attempts)
  pick <- function(field, type) {
    vapply(attempts, function(a) a[[field]] %||% NA, type)
  }
  out <- .sg_meta_frame(spec, ix)
  out$fit_id <- ix$fit_key
  out$status <- pick("status", character(1))
  out$error <- pick("error", character(1))
  out$elapsed <- pick("elapsed", numeric(1))
  out$warnings <- pick("warnings", character(1))
  rownames(out) <- NULL
  .sg_order_columns(.sg_frame(as.list(out), n))
}

# `measurements` holds, for each attempt, the measurement list of its record.
.sg_measure_frame <- function(spec, ix, measurements) {
  # Measurement order within each attempt; one pass, no growing vectors.
  declared <- names(spec$cap$measures)
  per <- lapply(measurements, function(m) {
    m <- m[declared[declared %in% names(m)]]
    lapply(m, `[[`, "rows")
  })
  counts <- lengths(per)
  frames <- unlist(per, recursive = FALSE, use.names = FALSE)
  which <- unlist(lapply(per, base::names), use.names = FALSE)
  attempt <- rep.int(seq_along(per), counts)
  n <- vapply(frames, nrow, integer(1))
  if (!sum(n)) {
    return(data.frame())
  }
  owner <- rep.int(attempt, n)
  which <- rep.int(which, n)
  values <- .sg_rows(frames[n > 0L])
  sub <- lapply(ix, function(x) x[owner])
  meta <- .sg_meta_frame(spec, sub)
  meta$measurement <- which
  out <- c(as.list(meta), as.list(values))
  .sg_order_columns(.sg_frame(out, length(owner)))
}

# Bind measurement rows. Rows of one measurement usually share their
# columns, so bind each column with one vec_c() when every frame has the same
# names; otherwise fall back to vec_rbind().
.sg_rows <- function(frames) {
  if (length(frames) == 1L) {
    return(frames[[1L]])
  }
  columns <- names(frames[[1L]])
  same <- all(vapply(frames, function(f) identical(names(f), columns), logical(1)))
  if (!same) {
    return(do.call(vctrs::vec_rbind, unname(frames)))
  }
  out <- lapply(columns, function(col) {
    vctrs::vec_c(!!!lapply(frames, `[[`, col))
  })
  names(out) <- columns
  .sg_frame(out, vctrs::vec_size(out[[1L]]))
}

# Comparisons ------------------------------------------------------------------------
#
# Group the attempts by each comparison's `by` columns and call its compute
# function once per group. Groups come from attempts, so a group whose fits
# all failed still appears. Returns a list of result frames by comparison.

.sg_compare <- function(
  comparisons,
  study_name,
  results,
  attempts,
  bundles = NULL,
  stream
) {
  attempts <- attempts[
    order(attempts$dataset_id, attempts$method_id, method = "radix"),
    ,
    drop = FALSE
  ]
  rownames(attempts) <- NULL
  if (nrow(results)) {
    results <- results[
      order(
        results$dataset_id,
        results$method_id,
        results$measurement,
        method = "radix"
      ),
      ,
      drop = FALSE
    ]
    rownames(results) <- NULL
    result_pair <- paste(results$dataset_id, results$method_id, sep = "/")
  }
  out <- list()
  for (name in names(comparisons)) {
    spec <- comparisons[[name]]
    grouping <- attempts[, spec$by, drop = FALSE]
    groups <- group_ids(grouping, spec$by)
    rows_out <- list()
    for (group in unique(groups)) {
      in_group <- groups == group
      members <- attempts[in_group, , drop = FALSE]
      rownames(members) <- NULL
      first <- which(in_group)[1L]
      rows <- if (nrow(results)) {
        results[
          result_pair %in% paste(members$dataset_id, members$method_id, sep = "/"),
          ,
          drop = FALSE
        ]
      } else {
        results
      }
      rownames(rows) <- NULL
      group_values <- if (length(spec$by)) {
        as.list(grouping[first, , drop = FALSE])
      } else {
        list()
      }
      artifacts <- list()
      if (!is.null(bundles)) {
        artifacts <- bundles[intersect(members$method_id, names(bundles))]
        for (id in names(artifacts)) {
          missing <- setdiff(spec$needs, names(artifacts[[id]]))
          if (length(missing)) {
            cli::cli_abort(c(
              "Comparison {.val {name}} needs discarded artifacts for method {.val {id}}: {paste(missing, collapse = ', ')}.",
              i = "Use a new store path to fit again with these artifacts kept."
            ))
          }
          artifacts[[id]] <- artifacts[[id]][spec$needs]
        }
      }
      context <- list(study = study_name, group = group_values, attempts = members)
      value <- .invoke(
        c(
          stream,
          name,
          paste(names(group_values), vapply(group_values, format, character(1)), sep = "=", collapse = "\r")
        ),
        spec$compute,
        results = rows,
        artifacts = artifacts,
        context = context
      )
      if (!is.data.frame(value)) {
        cli::cli_abort("Comparison {.val {name}} must return a data frame.")
      }
      if (nrow(value)) {
        if (any(names(value) %in% c("comparison", spec$by))) {
          cli::cli_abort("Comparison output uses reserved grouping columns.")
        }
        meta <- grouping[rep(first, nrow(value)), , drop = FALSE]
        meta$comparison <- rep(name, nrow(value))
        rows_out[[length(rows_out) + 1L]] <- .sg_frame(
          c(as.list(meta), as.list(value)),
          nrow(value)
        )
      }
    }
    out[[name]] <- .sg_bind(rows_out)
  }
  out
}
