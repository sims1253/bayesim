# Study jobs -------------------------------------------------------------------
#
# A job is one unit of work that a backend runs: prepare a condition, fit one
# method to one dataset (or measure an existing fit again), or run the
# artifact comparisons of one dataset. Workers read inputs from the store and
# write their records to it. They return only a short status and the name of
# the record they used, so the controller never receives fits or artifacts.
#
# The worker that runs a job holds its claim, from before it starts until
# after it has published, and only that worker releases it. The controller
# holds no claims. A worker that finds the claim held returns "busy"; one
# that finds the job done after it has claimed it returns the record and
# `ran = FALSE`.
#
# The worker holds the run specification serialised. A backend sends it to
# each process once; each job is then a short list of indices and keys.
# Each job unserialises its own copy, so that no job can see what another
# did to an environment that the study's functions capture (invariant 5).

.make_worker <- function(spec) {
  payload <- serialize(spec, NULL)
  rm(spec)
  cache <- new.env(parent = emptyenv())
  function(job) .run_job(job, unserialize(payload), cache)
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
  st <- .job_store(spec)
  token <- NULL
  on.exit(.release(st, job$key, token), add = TRUE)
  token <- .claim(st, job$key)
  if (is.null(token)) {
    return(list(kind = job$kind, status = "busy", ran = FALSE))
  }
  done <- .job_done(job, spec, st)
  if (!is.null(done)) {
    done$kind <- job$kind
    done$ran <- FALSE
    return(done)
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
  out$ran <- TRUE
  out$seconds <- proc.time()[["elapsed"]] - started
  out$warnings <- length(cache$warnings)
  out
}

# The output of a job that another process finished, or NULL. A failed
# attempt counts only for a new fit and only without `retry`.
.job_done <- function(job, spec, st) {
  if (job$kind == "prepare") {
    if (.store_has(.store_file(st, "prepared", job$key))) {
      return(list(status = "success"))
    }
  } else if (job$kind == "fit") {
    record <- paste0(job$fit_key, "-", job$mhash)
    if (.store_has(.store_file(st, "fits", record))) {
      return(list(status = "success", record = record))
    }
    if (job$mode == "fit" && !isTRUE(spec$retry)) {
      failures <- .fit_failures(st, job$fit_key)
      if (length(failures)) {
        return(list(status = "failed", record = .latest(failures)))
      }
    }
  } else if (all(.store_has(.store_file(st, "comparisons", job$keys)))) {
    return(list(status = "success"))
  }
  NULL
}

# Attempt IDs sort by time, so the last name is the latest attempt.
.latest <- function(names) {
  sort(names, method = "radix")[[length(names)]]
}

# Inputs ----------------------------------------------------------------------
#
# The caches of a worker hold serialised values. Every job, the first
# included, gets its own copy, so that a job that changes a dataset or a
# prepared value in place does not change what the next job sees.

.job_store <- function(spec) .store_handle(spec$root)

.job_prepared <- function(spec, cache, cidx) {
  if (is.null(spec$cap$prepare)) {
    return(NULL)
  }
  key <- spec$prepare_keys[[cidx]]
  cache$prepared <- cache$prepared %||% list()
  raw <- cache$prepared[[key]]
  if (is.null(raw)) {
    file <- .store_file(.job_store(spec), "prepared", key)
    if (!.store_has(file)) {
      cli::cli_abort(
        "The prepared value for condition {.val {spec$condition_ids[[cidx]]}} is missing from the store."
      )
    }
    raw <- .read_raw(file)
    # Keep the last few conditions; a worker usually sees one at a time.
    cache$prepared <- utils::tail(
      c(cache$prepared, stats::setNames(list(raw), key)),
      8L
    )
  }
  unserialize(raw)
}

.job_context <- function(job, spec, cache) {
  list(
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
# or generated again on its own stream. When the data are retained, a
# generated dataset is published at once, and a process that loses the race
# uses the dataset that won. Every fit on one dataset thus sees the stored
# dataset.
.job_data <- function(job, context, spec, cache) {
  if (!identical(cache$data_id, job$dataset_id)) {
    cache$data_id <- NULL
    st <- .job_store(spec)
    file <- .store_file(st, "datasets", job$dataset_id)
    raw <- if (.store_has(file)) {
      .read_raw(file)
    } else {
      cache$stage <- "generate"
      value <- .invoke(
        c("data", job$dataset_id),
        spec$cap$generate,
        condition = context$condition,
        context = context
      )
      raw <- serialize(value, NULL, xdr = FALSE)
      rm(value)
      if ("data" %in% spec$retention) {
        raw <- .publish_raw(raw, file, spec$compress)
      }
      raw
    }
    cache$data_raw <- raw
    cache$data_id <- job$dataset_id
  }
  unserialize(cache$data_raw)
}

# Prepare ----------------------------------------------------------------------

.job_prepare <- function(job, spec, cache) {
  condition <- spec$conditions[[job$cidx]]
  key <- spec$prepare_keys[[job$cidx]]
  value <- .invoke(
    c("prepare", key),
    spec$cap$prepare,
    condition = condition,
    context = list(condition = condition)
  )
  .write_payload(
    value,
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
  record <- paste0(job$fit_key, "-", job$mhash)

  if (identical(job$mode, "remeasure")) {
    return(.job_remeasure(job, spec, cache, st, method, context, get_data, record))
  }

  attempt_id <- .attempt_id()
  get_data()
  fit <- .job_run_fit(job, spec, cache, method, context, data[[1L]])
  attempt <- list(
    id = attempt_id,
    status = if (is.null(fit$error)) "success" else "failed",
    error = fit$error %||% NA_character_,
    elapsed = fit$elapsed,
    warnings = .job_warnings(cache)
  )
  if (!is.null(fit$error)) {
    failed <- .failure_record(job$fit_key, attempt_id)
    .write_payload(
      list(attempt = attempt, extract_id = method$extract_id),
      .store_file(st, "fits", failed)
    )
    return(list(status = "failed", record = failed))
  }
  artifacts <- .job_extract(
    job, spec, cache, method, context, fit$value, data[[1L]],
    spec$required[spec$required != "data"]
  )
  rm(fit)
  measurements <- .job_measure(
    job, spec, cache, context, list(),
    function(names) {
      if ("data" %in% names) artifacts["data"] <<- list(get_data())
      artifacts[names]
    }
  )
  # The data are a dataset, not an artifact of the attempt.
  kept <- setdiff(
    intersect(names(artifacts), c(spec$retention, spec$comparison_needs)),
    "data"
  )
  files <- .job_publish_artifacts(st, spec, job$fit_key, attempt_id, artifacts[kept])$files
  attempt$warnings <- .job_warnings(cache, attempt$warnings)
  won <- .write_payload(
    list(
      attempt = attempt,
      measurements = measurements,
      extract_id = method$extract_id
    ),
    .store_file(st, "fits", record)
  )
  if (!won) {
    # Another process published this record first. No record refers to the
    # artifacts of this attempt.
    .store_delete(files)
  }
  list(status = "success", record = record)
}

# Measure an existing attempt again. The run planned the success records of
# this fit; their canonical attempt is the one with the smallest ID. Every
# measurement of that attempt is reused, whichever record holds it. An
# artifact that is not in the store is computed again: the fit runs again on
# its own stream, which gives the same fit when the fit depends only on its
# seed.
.job_remeasure <- function(job, spec, cache, st, method, context, get_data, record) {
  olds <- lapply(job$old, function(r) .read_payload(.store_file(st, "fits", r)))
  olds <- Filter(function(o) identical(o$extract_id, method$extract_id), olds)
  if (!length(olds)) {
    cli::cli_abort(c(
      "Extractor recipe changed for cached method {.val {context$method_id}}.",
      i = "Change the method version to fit again."
    ))
  }
  ids <- vapply(olds, function(o) o$attempt$id %||% "", character(1))
  canonical <- sort(unique(ids), method = "radix")[[1L]]
  olds <- olds[ids == canonical]
  attempt <- olds[[1L]]$attempt
  cached <- list()
  for (o in olds) {
    for (name in names(o$measurements)) {
      current <- spec$cap$measures[[name]]
      if (!is.null(current) && identical(o$measurements[[name]]$id, current$id)) {
        cached[[name]] <- o$measurements[[name]]
      }
    }
  }
  keep <- c(spec$retention, spec$comparison_needs)
  artifacts <- list()
  # Load or compute again the artifacts `names`.
  supply <- function(names) {
    names <- setdiff(names, names(artifacts))
    if ("data" %in% names) {
      artifacts["data"] <<- list(get_data())
    }
    names <- setdiff(names, "data")
    if (length(names)) {
      got <- .job_stored_artifacts(st, spec, job$fit_key, canonical, names)
      missing <- setdiff(names, names(got))
      if (length(missing)) {
        fresh <- .job_regenerate(job, spec, cache, method, context, get_data(), missing)
        publish <- intersect(missing, keep)
        fresh[publish] <- .job_publish_artifacts(
          st, spec, job$fit_key, canonical, fresh[publish]
        )$artifacts
        got[missing] <- fresh[missing]
      }
      artifacts[names(got)] <<- got
    }
    artifacts
  }
  measurements <- .job_measure(
    job, spec, cache, context, cached,
    function(names) supply(names)[names]
  )
  # The record promises the retained artifacts and data. Put them in the
  # store first.
  retained <- setdiff(spec$retention, names(artifacts))
  stored <- .store_has(vapply(
    setdiff(retained, "data"),
    function(n) .artifact_file(st, "artifacts", job$fit_key, canonical, n),
    character(1)
  ))
  absent <- setdiff(retained, "data")[!stored]
  if (length(absent)) {
    supply(absent)
  }
  if ("data" %in% retained) {
    supply("data")
  }
  attempt$warnings <- .job_warnings(cache, attempt$warnings)
  .write_payload(
    list(
      attempt = attempt,
      measurements = measurements,
      extract_id = method$extract_id
    ),
    .store_file(st, "fits", record)
  )
  list(status = "success", record = record)
}

# Run the method's fit on `data`. Returns list(value, elapsed) or
# list(error, elapsed).
.job_run_fit <- function(job, spec, cache, method, context, data) {
  cache$stage <- "fit"
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    list(
      value = .invoke(
        c("fit", job$dataset_id, context$method_id),
        method$fit,
        data = data,
        context = context
      )
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  fit$elapsed <- proc.time()[["elapsed"]] - started
  fit
}

.job_extract <- function(job, spec, cache, method, context, fit, data, names) {
  cache$stage <- "extract"
  artifacts <- list()
  for (name in names) {
    artifacts[name] <- list(
      if (name == "fit") {
        fit
      } else {
        .invoke(
          c("extract", job$dataset_id, context$method_id, name),
          method$extract[[name]],
          fit = fit,
          data = data,
          context = context
        )
      }
    )
  }
  artifacts
}

# Compute again the artifacts `names` of a fit that succeeded before.
.job_regenerate <- function(job, spec, cache, method, context, data, names) {
  fit <- .job_run_fit(job, spec, cache, method, context, data)
  if (!is.null(fit$error)) {
    cli::cli_abort(c(
      "The fit failed when it ran again to compute discarded artifacts: {fit$error}",
      i = "The fit succeeded before on the same seed, so it does not depend on its seed alone."
    ))
  }
  .job_extract(job, spec, cache, method, context, fit$value, data, names)
}

# The artifacts `names` of one attempt that the store holds. A retained
# artifact that only the transient files hold is kept before it is used, so
# that a record never names a kept artifact that is not there.
.job_stored_artifacts <- function(st, spec, fit_key, attempt_id, names) {
  got <- .job_load_artifacts(st, fit_key, attempt_id, names)
  promote <- intersect(attr(got, "transient"), spec$retention)
  attr(got, "transient") <- NULL
  if (length(promote)) {
    got[promote] <- .job_publish_artifacts(
      st, spec, fit_key, attempt_id, got[promote]
    )$artifacts
  }
  got
}

# The artifacts `names` of one attempt that the store holds, with the names
# found in the transient files in the attribute "transient".
.job_load_artifacts <- function(st, fit_key, attempt_id, names) {
  out <- list()
  transient <- character()
  for (name in names) {
    for (kind in c("artifacts", "transient")) {
      file <- .artifact_file(st, kind, fit_key, attempt_id, name)
      if (!.store_has(file)) {
        next
      }
      # A comparison job may delete a transient artifact meanwhile.
      raw <- tryCatch(.read_raw(file), error = function(e) {
        if (file.exists(file)) stop(e)
        NULL
      })
      if (!is.null(raw)) {
        out[name] <- list(unserialize(raw))
        if (kind == "transient") {
          transient <- c(transient, name)
        }
        break
      }
    }
  }
  attr(out, "transient") <- transient
  out
}

# Publish artifacts of one attempt: retained ones to `artifacts`, the ones
# that only comparisons need to `transient`. When another process published
# a retained artifact of the same attempt first, its value replaces this
# one, so that the measurements agree with the kept file. A comparison job
# may delete a transient file at any time, so a process that loses the race
# to publish one keeps its own value. Losing a race is never an error.
.job_publish_artifacts <- function(st, spec, fit_key, attempt_id, artifacts) {
  files <- character()
  for (name in names(artifacts)) {
    kind <- if (name %in% spec$retention) "artifacts" else "transient"
    file <- .artifact_file(st, kind, fit_key, attempt_id, name)
    raw <- serialize(artifacts[[name]], NULL, xdr = FALSE)
    if (.write_raw(raw, file, spec$compress)) {
      files <- c(files, file)
    } else if (kind == "artifacts") {
      artifacts[name] <- list(.read_payload(file))
    }
  }
  list(artifacts = artifacts, files = files)
}

# Compute the measurements that `cached` does not hold. `supply(names)`
# returns the artifacts `names`.
.job_measure <- function(job, spec, cache, context, cached, supply) {
  measurements <- list()
  for (name in names(spec$cap$measures)) {
    measure <- spec$cap$measures[[name]]
    old <- cached[[name]]
    if (!is.null(old) && identical(old$id, measure$id)) {
      measurements[[name]] <- old
      next
    }
    artifacts <- supply(measure$needs)
    cache$stage <- "measure"
    rows <- .invoke(
      c("measure", job$dataset_id, context$method_id, name),
      measure$compute,
      artifacts = artifacts,
      context = context
    )
    measurements[[name]] <- list(
      id = measure$id,
      rows = .job_check_rows(rows, name, spec)
    )
  }
  measurements
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
# fitted, store each result, then delete the transient artifacts. The
# controller passes the records it uses and their attempt IDs; the artifacts
# are those of these attempts, computed again when they are gone.
.job_compare <- function(job, spec, cache) {
  st <- .job_store(spec)
  records <- lapply(job$records, function(r) {
    .read_payload(.store_file(st, "fits", r))
  })
  ids <- vapply(records, function(r) r$attempt$id %||% "", character(1))
  if (!identical(ids, job$attempts)) {
    cli::cli_abort("A record of dataset {.val {job$dataset_id}} changed while the run used it.")
  }
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
  needs <- setdiff(spec$comparison_needs, "data")
  bundles <- list()
  for (i in which(status == "success")) {
    midx <- job$midx[[i]]
    a <- .job_stored_artifacts(st, spec, job$fit_keys[[i]], ids[[i]], needs)
    missing <- setdiff(needs, names(a))
    if (length(missing)) {
      method <- spec$cap$methods[[midx]]
      mctx <- .job_method_context(context, spec, midx)
      fresh <- .job_regenerate(
        job, spec, cache, method, mctx,
        .job_data(job, context, spec, cache), missing
      )
      # Retained artifacts go back to the store; the others are not kept.
      publish <- intersect(missing, spec$retention)
      fresh[publish] <- .job_publish_artifacts(
        st, spec, job$fit_keys[[i]], ids[[i]], fresh[publish]
      )$artifacts
      a[missing] <- fresh[missing]
    }
    if ("data" %in% spec$comparison_needs) {
      a["data"] <- list(.job_data(job, context, spec, cache))
    }
    bundles[[spec$method_ids[[midx]]]] <- a
  }
  cache$stage <- "compare"
  for (name in names(job$keys)) {
    file <- .store_file(st, "comparisons", job$keys[[name]])
    if (.store_has(file)) {
      next
    }
    value <- .sg_compare(
      spec$cap$comparisons[name],
      results,
      attempts,
      bundles = bundles,
      stream = c("compare", job$dataset_id)
    )
    .write_payload(list(value = value[[name]]), file)
  }
  # A retained artifact that was transient is now kept, so every transient
  # file of these attempts can go.
  transient <- needs
  if (length(transient)) {
    .store_delete(unlist(lapply(which(status == "success"), function(i) {
      vapply(transient, function(name) {
        .artifact_file(st, "transient", job$fit_keys[[i]], ids[[i]], name)
      }, character(1))
    })))
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
      context <- list(group = group_values, attempts = members)
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
