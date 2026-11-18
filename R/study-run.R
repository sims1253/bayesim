# Study runs -----------------------------------------------------------------------
#
# run_study() derives every key from the declaration, lists the store, and
# runs only the jobs whose records are missing. There is no manifest: the
# plan is the set difference between the keys the study needs and the file
# names in the store.

#' Run a study
#'
#' `run_study()` generates the datasets, fits every method, computes the
#' measurements and comparisons, and returns the results. Each fit is written
#' to the store as soon as it finishes. A later call with the same store runs
#' only the work that is missing: new replicates, new methods, or new
#' measurements on kept artifacts.
#'
#' Results do not depend on the backend, the number of workers or the order
#' in which jobs finish. Several R processes can run the same study on one
#' store at the same time; each job runs once.
#'
#' An error in a method's `fit` is recorded as a failed attempt and the study
#' continues. A failed attempt is not tried again unless `retry = TRUE`. An
#' error in any other function stops the run. Records that finished before the
#' error stay in the store.
#'
#' @inheritParams plan_study
#' @param path Store directory. With `NULL`, the run uses a temporary store
#'   that is deleted when the run ends.
#' @param backend Where jobs run: [sequential_backend()] (the default) or
#'   [mirai_backend()].
#' @param retry If `TRUE`, fit failed attempts again.
#' @param progress If `TRUE`, show a progress bar.
#'
#' @return A `bayesim_study_run`: a list with the data frames `attempts`,
#'   `measurements` and `comparisons`, the `study`, the `plan`, the `path`
#'   and `executed`, a data frame of the jobs that this call ran.
#' @seealso [summarise_study()], [assess_study()]
#' @export
run_study <- function(
  study,
  replicates,
  seed = 1L,
  path = NULL,
  backend = sequential_backend(),
  retry = FALSE,
  progress = interactive()
) {
  if (!inherits(backend, "bayesim_backend")) {
    cli::cli_abort(
      "{.arg backend} must be {.fn sequential_backend} or {.fn mirai_backend}."
    )
  }
  setup <- .sg_setup(study, replicates, seed)
  if (is.null(path)) {
    root <- tempfile("bayesim-study-")
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
  } else {
    root <- path
  }
  st <- .store_open(root)
  setup$spec$root <- st$root
  executed <- .sg_execute(setup, st, backend, retry, isTRUE(progress))
  run <- .sg_collect(setup, st, seq_len(nrow(setup$datasets)), study)
  run$path <- path
  run$executed <- executed
  run
}

#' Run one dataset of a study
#'
#' `evaluate_replicate()` runs all methods on the dataset of one condition and
#' replicate, in this R process. An external scheduler can call it for each
#' dataset, with one shared store. The records are the same that
#' [run_study()] writes, so either function can continue the other's work.
#'
#' @inheritParams run_study
#' @param condition_id A condition ID of the study.
#' @param replicate A replicate number.
#'
#' @return A list with `attempts`, `measurements` and `comparisons`. The
#'   comparisons are the artifact comparisons only; comparisons of saved
#'   measurements need all datasets, so use [assess_study()] on a full run
#'   for them.
#' @export
evaluate_replicate <- function(
  study,
  condition_id,
  replicate,
  seed = 1L,
  path = NULL
) {
  .sg_study(study)
  .sg_count(replicate, "replicate")
  cidx <- match(condition_id, study$conditions$condition_id)
  if (length(condition_id) != 1L || is.na(cidx)) {
    cli::cli_abort("Unknown {.arg condition_id}: {.val {condition_id}}.")
  }
  setup <- .sg_setup(study, replicate, seed)
  d <- which(
    setup$datasets$cidx == cidx & setup$datasets$replicate == replicate
  )
  setup <- .sg_restrict(setup, d)
  if (is.null(path)) {
    root <- tempfile("bayesim-study-")
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
  } else {
    root <- path
  }
  st <- .store_open(root)
  setup$spec$root <- st$root
  .sg_execute(setup, st, sequential_backend(), retry = FALSE, progress = FALSE)
  run <- .sg_collect(setup, st, 1L, study, saved = FALSE)
  run[c("attempts", "measurements", "comparisons")]
}

#' Apply an analysis policy to a study run
#'
#' A policy keeps or drops measurement rows, for example fits with poor
#' convergence. `assess_study()` then runs the comparisons of saved
#' measurements again on the rows that remain. Attempts, including failed
#' ones, stay visible.
#'
#' @param run A [run_study()] result.
#' @param policy `NULL` to keep every row, or
#'   `function(measurements, attempts)` that returns one `TRUE` or `FALSE`
#'   per measurement row.
#'
#' @return A `bayesim_study_assessment`: a list with `measurements` (kept
#'   rows), `excluded`, `attempts`, `comparisons` (computed on the kept rows),
#'   `execution_comparisons` (the stored artifact comparisons) and `study`.
#' @export
assess_study <- function(run, policy = NULL) {
  if (!inherits(run, "bayesim_study_run")) {
    cli::cli_abort("{.arg run} must be a {.fn run_study} result.")
  }
  rows <- run$measurements
  keep <- if (is.null(policy)) rep(TRUE, nrow(rows)) else policy(rows, run$attempts)
  if (!is.logical(keep) || length(keep) != nrow(rows) || anyNA(keep)) {
    cli::cli_abort(
      "{.arg policy} must return one nonmissing logical value per measurement row."
    )
  }
  selected <- rows[keep, , drop = FALSE]
  saved <- Filter(function(x) !length(x$needs), run$study$comparisons)
  raw_names <- setdiff(names(run$study$comparisons), names(saved))
  structure(
    list(
      measurements = selected,
      excluded = rows[!keep, , drop = FALSE],
      attempts = run$attempts,
      comparisons = .sg_bind(.sg_compare(
        saved,
        selected,
        run$attempts,
        stream = c("compare", run$plan$seed)
      )),
      execution_comparisons = if (nrow(run$comparisons)) {
        run$comparisons[
          run$comparisons$comparison %in% raw_names,
          ,
          drop = FALSE
        ]
      } else {
        run$comparisons
      },
      study = run$study
    ),
    class = "bayesim_study_assessment"
  )
}

# Setup --------------------------------------------------------------------------------

# Derive every key of the study. Datasets run by replicate, then condition;
# fits by dataset, then method in study order.
.sg_setup <- function(study, replicates, seed) {
  plan <- .sg_plan(study, replicates, seed)
  cap <- attr(plan, "capture")
  attr(plan, "capture") <- NULL
  seed <- as.integer(seed)
  replicates <- as.integer(replicates)
  conditions <- study$conditions
  n_cond <- nrow(conditions)
  n_method <- length(cap$methods)
  method_ids <- names(cap$methods)
  rows <- lapply(seq_len(n_cond), function(i) as.list(conditions[i, , drop = FALSE]))
  cond_keys <- .condition_keys(rows)
  prepare_keys <- if (is.null(cap$prepare)) {
    rep("", n_cond)
  } else {
    .hash_text("prepare", seed, cond_keys, cap$prepare_id)
  }
  cidx <- rep(seq_len(n_cond), replicates)
  replicate <- rep(seq_len(replicates), each = n_cond)
  datasets <- .sg_frame(
    list(
      cidx = cidx,
      replicate = replicate,
      dataset_id = .hash_text(
        "dataset",
        seed,
        cond_keys[cidx],
        replicate,
        prepare_keys[cidx],
        cap$generate_id
      )
    ),
    length(cidx)
  )
  measure_text <- paste(
    sort(paste(names(cap$measures), vapply(cap$measures, `[[`, "", "id"), sep = ":")),
    collapse = ","
  )
  fit_ids <- vapply(cap$methods, `[[`, "", "fit_id")
  mhash <- .hash_text(
    "measures",
    vapply(cap$methods, `[[`, "", "extract_id"),
    measure_text,
    paste(sort(cap$retention), collapse = ",")
  )
  didx <- rep(seq_len(nrow(datasets)), each = n_method)
  midx <- rep(seq_len(n_method), nrow(datasets))
  fits <- .sg_frame(
    list(
      didx = didx,
      midx = midx,
      fit_key = .hash_text(
        "fit",
        datasets$dataset_id[didx],
        method_ids[midx],
        fit_ids[midx]
      )
    ),
    length(didx)
  )
  fits$record <- paste0(fits$fit_key, "-", mhash[midx])
  raw <- Filter(function(x) length(x$needs) > 0L, cap$comparisons)
  cond_meta <- .sg_condition_meta(conditions)
  method_meta <- .sg_method_meta(study$methods)
  spec <- list(
    cap = cap,
    root = NULL,
    seed = seed,
    conditions = rows,
    condition_ids = conditions$condition_id,
    method_ids = method_ids,
    prepare_keys = prepare_keys,
    cond_meta = cond_meta,
    method_meta = method_meta,
    meta_names = c(
      "condition_id",
      "dataset_id",
      "replicate",
      "method_id",
      names(cond_meta),
      names(method_meta)
    ),
    required = cap$required,
    retention = cap$retention,
    comparison_needs = unique(unlist(lapply(raw, `[[`, "needs"), use.names = FALSE)),
    raw_comparisons = names(raw),
    raw_ids = vapply(raw, `[[`, "", "id"),
    compress = cap$compress
  )
  list(plan = plan, spec = spec, datasets = datasets, fits = fits, mhash = mhash)
}

.sg_restrict <- function(setup, d) {
  keep <- setup$fits$didx %in% d
  setup$fits <- setup$fits[keep, , drop = FALSE]
  setup$fits$didx <- match(setup$fits$didx, d)
  setup$datasets <- setup$datasets[d, , drop = FALSE]
  rownames(setup$fits) <- NULL
  rownames(setup$datasets) <- NULL
  setup
}

# The state of each fit, from the store's file names: "done" (a record with
# the current measurements), "failed", "measure" (a successful record with
# other measurements) or "fit".
.sg_fit_states <- function(setup, names, retry = FALSE) {
  fits <- setup$fits
  ok <- names[!endsWith(names, ".failed")]
  ok_keys <- substr(ok, 1L, 40L)
  failed_keys <- substr(names[endsWith(names, ".failed")], 1L, 40L)
  state <- rep("fit", nrow(fits))
  hit <- match(fits$fit_key, ok_keys)
  state[!is.na(hit)] <- "measure"
  if (!retry) {
    state[state == "fit" & fits$fit_key %in% failed_keys] <- "failed"
  }
  state[fits$record %in% ok] <- "done"
  fits$state <- state
  fits$old <- ok[hit]
  fits
}

# The key of each artifact comparison of one dataset. It depends on the
# records of every method, so a new fit, a retried failure or a new
# measurement gives a new key.
.sg_compare_keys <- function(setup, dataset_id, records) {
  spec <- setup$spec
  members <- paste(sort(records), collapse = ",")
  keys <- .hash_text("comparison", spec$raw_ids, dataset_id, members)
  names(keys) <- spec$raw_comparisons
  keys
}

# Execution -----------------------------------------------------------------------------

.sg_execute <- function(setup, st, backend, retry, progress) {
  # Backends may draw random numbers, for example mirai for its socket
  # names. The caller's stream must not move.
  withr::local_preserve_seed()
  spec <- setup$spec
  datasets <- setup$datasets
  index <- .store_index(st, c("fits", "prepared", "comparisons"))
  fits <- .sg_fit_states(setup, index$fits, retry = retry)
  n_fit <- nrow(fits)
  status <- ifelse(fits$state == "failed", "failed", "success")
  status[fits$state %in% c("fit", "measure")] <- NA
  records <- fits$record
  records[fits$state == "failed"] <- paste0(fits$fit_key[fits$state == "failed"], ".failed")
  todo <- which(fits$state %in% c("fit", "measure"))
  pending <- tabulate(fits$didx[todo], nrow(datasets))
  has_raw <- length(spec$raw_comparisons) > 0L
  done_comparisons <- index$comparisons

  log <- list(kind = character(), id = character())
  n_log <- 0L
  held <- character()
  handle <- NULL
  on.exit({
    if (!is.null(handle)) handle$close()
    for (key in held) .release(st, key)
  }, add = TRUE)

  # Queues: `front` holds prepare and compare jobs, which run first; `queue`
  # holds fit indices in study order.
  front <- list()
  queue <- todo
  qpos <- 1L
  conditions_todo <- unique(datasets$cidx[fits$didx[todo]])
  blocked <- rep(FALSE, length(spec$condition_ids))
  if (!is.null(spec$cap$prepare)) {
    for (c in conditions_todo) {
      if (!spec$prepare_keys[[c]] %in% index$prepared) {
        blocked[c] <- TRUE
        front[[length(front) + 1L]] <- list(
          kind = "prepare",
          cidx = c,
          key = spec$prepare_keys[[c]]
        )
      }
    }
  }
  waiting <- vector("list", length(blocked))

  compare_job <- function(d) {
    rows <- which(fits$didx == d)
    keys <- .sg_compare_keys(setup, datasets$dataset_id[[d]], records[rows])
    keys <- keys[!keys %in% done_comparisons]
    if (!length(keys)) {
      return(NULL)
    }
    list(
      kind = "compare",
      key = paste0("cmp-", datasets$dataset_id[[d]]),
      d = d,
      cidx = datasets$cidx[[d]],
      replicate = datasets$replicate[[d]],
      dataset_id = datasets$dataset_id[[d]],
      midx = fits$midx[rows],
      fit_keys = fits$fit_key[rows],
      records = records[rows],
      keys = keys
    )
  }
  if (has_raw) {
    for (d in which(pending == 0L)) {
      job <- compare_job(d)
      if (!is.null(job)) front[[length(front) + 1L]] <- job
    }
  }

  fit_job <- function(i) {
    d <- fits$didx[[i]]
    list(
      kind = "fit",
      key = fits$fit_key[[i]],
      i = i,
      mode = if (fits$state[[i]] == "measure") "remeasure" else "fit",
      old = fits$old[[i]],
      cidx = datasets$cidx[[d]],
      replicate = datasets$replicate[[d]],
      dataset_id = datasets$dataset_id[[d]],
      midx = fits$midx[[i]],
      fit_key = fits$fit_key[[i]],
      mhash = setup$mhash[[fits$midx[[i]]]]
    )
  }

  next_job <- function() {
    if (length(front)) {
      job <- front[[1L]]
      front[[1L]] <<- NULL
      return(job)
    }
    while (qpos <= length(queue)) {
      i <- queue[[qpos]]
      qpos <<- qpos + 1L
      c <- datasets$cidx[[fits$didx[[i]]]]
      if (blocked[[c]]) {
        waiting[[c]] <<- c(waiting[[c]], i)
        next
      }
      return(fit_job(i))
    }
    NULL
  }

  # What a finished job, or one that another process finished, changes.
  complete <- function(job, result) {
    if (job$kind == "prepare") {
      blocked[[job$cidx]] <<- FALSE
      if (length(waiting[[job$cidx]])) {
        queue <<- c(waiting[[job$cidx]], queue[seq.int(qpos, length.out = length(queue) - qpos + 1L)])
        qpos <<- 1L
        waiting[job$cidx] <<- list(NULL)
      }
    } else if (job$kind == "fit") {
      status[[job$i]] <<- result$status
      records[[job$i]] <<- result$record
      d <- fits$didx[[job$i]]
      pending[[d]] <<- pending[[d]] - 1L
      if (has_raw && pending[[d]] == 0L) {
        cj <- compare_job(d)
        if (!is.null(cj)) front <<- c(list(cj), front)
      }
    } else {
      done_comparisons <<- c(done_comparisons, job$keys)
    }
    if (!is.null(pb)) cli::cli_progress_update(id = pb)
  }

  # The output of a job, if it exists: list(status, record) or NULL.
  finished <- function(job) {
    if (job$kind == "prepare") {
      if (.store_has(.store_file(st, "prepared", job$key))) {
        return(list(status = "success"))
      }
    } else if (job$kind == "fit") {
      record <- fits$record[[job$i]]
      if (.store_has(.store_file(st, "fits", record))) {
        return(list(status = "success", record = record))
      }
      failed <- paste0(job$fit_key, ".failed")
      if (!retry && .store_has(.store_file(st, "fits", failed))) {
        return(list(status = "failed", record = failed))
      }
    } else if (all(.store_has(.store_file(st, "comparisons", job$keys)))) {
      return(list(status = "success"))
    }
    NULL
  }

  total <- length(front) + length(todo)
  pb <- NULL
  if (progress && total) {
    pb <- cli::cli_progress_bar("Running study", total = total, clear = FALSE)
  }

  handle <- backend$open(.make_worker(spec))
  capacity <- handle$capacity
  inflight <- list()
  deferred <- list()
  warned <- 0L

  repeat {
    while (length(inflight) < capacity) {
      job <- next_job()
      if (is.null(job)) {
        break
      }
      if (!.claim(st, job$key)) {
        deferred[[length(deferred) + 1L]] <- job
        next
      }
      held <- c(held, job$key)
      done <- finished(job)
      if (!is.null(done)) {
        .release(st, job$key)
        held <- held[held != job$key]
        complete(job, done)
        next
      }
      n_log <- n_log + 1L
      id <- as.character(n_log)
      log$kind[[n_log]] <- job$kind
      log$id[[n_log]] <- job$key
      inflight[[id]] <- job
      handle$submit(id, job[!names(job) %in% c("key", "i", "d")])
    }
    if (!length(inflight)) {
      if (!length(deferred)) {
        break
      }
      # Another process holds the remaining jobs. Wait for it.
      Sys.sleep(0.05)
      still <- list()
      for (job in deferred) {
        if (dir.exists(file.path(st$claims, job$key)) &&
            !.claim_stale(file.path(st$claims, job$key))) {
          still[[length(still) + 1L]] <- job
          next
        }
        done <- finished(job)
        if (is.null(done)) {
          front <- c(front, list(job))
        } else {
          complete(job, done)
        }
      }
      deferred <- still
      next
    }
    for (res in handle$poll()) {
      job <- inflight[[res$id]]
      inflight[[res$id]] <- NULL
      value <- res$value
      if (!is.list(value) || is.null(value$kind)) {
        value <- list(error = list(
          message = paste("The worker failed:", paste(format(value), collapse = " ")),
          stage = job$kind,
          condition_id = spec$condition_ids[[job$cidx]],
          replicate = job$replicate,
          method_id = if (!is.null(job$midx) && job$kind == "fit") spec$method_ids[[job$midx]]
        ))
      }
      if (!is.null(value$error)) {
        .sg_stage_error(value$error)
      }
      .release(st, job$key)
      held <- held[held != job$key]
      complete(job, value)
      if (isTRUE(value$warnings > 0L)) {
        warned <- warned + 1L
      }
    }
  }
  if (!is.null(pb)) cli::cli_progress_done(id = pb)
  if (warned > 0L) {
    cli::cli_warn(
      "{warned} job{?s} raised warnings. See the {.field warnings} column of {.code attempts}."
    )
  }
  .sg_frame(log, length(log$kind))
}

.sg_stage_error <- function(e) {
  where <- c(
    if (!is.null(e$condition_id)) paste0("condition ", e$condition_id),
    if (!is.null(e$replicate)) paste0("replicate ", e$replicate),
    if (!is.null(e$method_id)) paste0("method ", e$method_id)
  )
  msg <- e$message
  cli::cli_abort(
    c("{msg}", i = "Stage {.field {e$stage}}{if (length(where)) paste0(', ', paste(where, collapse = ', ')) else ''}."),
    class = "bayesim_study_error",
    call = NULL
  )
}

# Collection -------------------------------------------------------------------------------

.sg_collect <- function(setup, st, d, study, saved = TRUE) {
  spec <- setup$spec
  datasets <- setup$datasets
  fits <- setup$fits[setup$fits$didx %in% d, , drop = FALSE]
  names <- .store_index(st, "fits")$fits
  record <- fits$record
  missing <- !record %in% names
  record[missing] <- paste0(fits$fit_key[missing], ".failed")
  files <- .store_file(st, "fits", record)
  sizes <- file.size(files)
  payloads <- lapply(seq_along(files), function(i) .read_payload(files[[i]], sizes[[i]]))
  ix <- list(
    dataset_id = datasets$dataset_id[fits$didx],
    cidx = datasets$cidx[fits$didx],
    replicate = datasets$replicate[fits$didx],
    midx = fits$midx,
    fit_key = fits$fit_key
  )
  attempts <- .sg_attempt_frame(spec, ix, lapply(payloads, `[[`, "attempt"))
  measurements <- .sg_measure_frame(
    spec,
    ix,
    lapply(payloads, function(p) p$measurements %||% list())
  )
  if (!nrow(measurements)) {
    measurements <- data.frame()
  }
  comparisons <- list()
  if (length(spec$raw_comparisons)) {
    for (dd in d) {
      rows <- which(fits$didx == dd)
      keys <- .sg_compare_keys(setup, datasets$dataset_id[[dd]], record[rows])
      for (k in keys) {
        comparisons[[length(comparisons) + 1L]] <- .read_payload(
          .store_file(st, "comparisons", k)
        )$value
      }
    }
  }
  if (saved) {
    plain <- Filter(function(x) !length(x$needs), study$comparisons)
    comparisons <- c(
      comparisons,
      .sg_compare(
        plain,
        measurements,
        attempts,
        stream = c("compare", spec$seed)
      )
    )
  }
  structure(
    list(
      study = study,
      plan = setup$plan,
      path = NULL,
      attempts = attempts,
      measurements = measurements,
      comparisons = .sg_bind(comparisons),
      executed = NULL
    ),
    class = "bayesim_study_run"
  )
}

#' @export
print.bayesim_study_run <- function(x, ...) {
  failed <- sum(x$attempts$status == "failed")
  cat(
    "<study run> ", x$study$name, ": ",
    nrow(x$attempts), " attempts (", failed, " failed), ",
    nrow(x$measurements), " measurement rows, ",
    nrow(x$comparisons), " comparison rows\n",
    sep = ""
  )
  invisible(x)
}
