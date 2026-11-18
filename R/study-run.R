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
#' store at the same time. A process claims each job before it runs it, so
#' that the others do other work, and usually each job runs once. A process
#' that finds a job claimed uses the other process's record as soon as it is
#' written. When a claim is older than `getOption("bayesim.claim_timeout")`
#' seconds (default one day), another process takes the job over. The
#' process that calls `run_study()` reads the option and passes it to its
#' workers. Taking over a job that is
#' still running costs only duplicate work: the first record written wins,
#' and every run uses it.
#'
#' An error in a method's `fit` is recorded as a failed attempt and the study
#' continues. A failed attempt is fitted again only by a run with
#' `retry = TRUE`; a run without it uses the latest failed attempt. A
#' successful attempt replaces failed ones for every later run. A run decides
#' at its start which record each fit uses, and reports those records, also
#' when another process retries a failed fit meanwhile. An error in any other
#' function stops the run. Records that finished before the error stay in the
#' store.
#'
#' A new measurement or comparison of an existing fit needs its artifacts.
#' When the store did not keep them, the fit runs again on its own stream;
#' this gives the same fit when the fit depends only on its seed.
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
  exec <- .sg_execute(setup, st, backend, retry, isTRUE(progress))
  run <- .sg_collect(setup, st, exec, study)
  run$path <- path
  run$executed <- exec$executed
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
  exec <- .sg_execute(setup, st, sequential_backend(), retry = FALSE, progress = FALSE)
  run <- .sg_collect(setup, st, exec, study, saved = FALSE)
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
# the current measurements), "failed" (only failed attempts, and no
# `retry`), "measure" (a successful record with other measurements) or
# "fit". `used` is the record that a done or failed fit uses: the current
# record, or the latest failed attempt. `old` lists, for each fit, its
# successful records.
.sg_fit_states <- function(setup, names, retry = FALSE) {
  fits <- setup$fits
  keys <- .fit_key_of(names)
  ok <- .is_success_record(names)
  ok_names <- names[ok]
  ok_keys <- keys[ok]
  failed <- sort(names[.is_failure_record(names)], decreasing = TRUE, method = "radix")
  latest_failure <- failed[match(fits$fit_key, .fit_key_of(failed))]
  state <- rep("fit", nrow(fits))
  state[fits$fit_key %in% ok_keys] <- "measure"
  if (!retry) {
    state[state == "fit" & !is.na(latest_failure)] <- "failed"
  }
  state[fits$record %in% ok_names] <- "done"
  fits$state <- state
  used <- rep(NA_character_, nrow(fits))
  used[state == "done"] <- fits$record[state == "done"]
  used[state == "failed"] <- latest_failure[state == "failed"]
  fits$used <- used
  old <- vector("list", nrow(fits))
  remeasure <- which(state == "measure")
  if (length(remeasure)) {
    by_key <- split(ok_names, ok_keys)
    old[remeasure] <- by_key[fits$fit_key[remeasure]]
  }
  fits$old <- old
  fits
}

# The members of a comparison: each record it reads with the attempt in it.
.sg_members <- function(records, attempts) {
  paste(sort(paste(records, attempts, sep = "@"), method = "radix"), collapse = ",")
}

# The key of each artifact comparison of one dataset. It names the records of
# every method and their attempts, so a new fit, a retried failure or a new
# measurement gives a new key. The name is in the key because it is in the
# comparison's seed and in its rows.
.sg_compare_keys <- function(setup, dataset_id, records, attempts) {
  spec <- setup$spec
  keys <- .hash_text(
    "comparison",
    spec$raw_comparisons,
    spec$raw_ids,
    dataset_id,
    .sg_members(records, attempts)
  )
  names(keys) <- spec$raw_comparisons
  keys
}

# The key of each comparison of saved measurements. It reads every record of
# the run.
.sg_saved_keys <- function(setup, comparisons, records, attempts) {
  keys <- .hash_text(
    "saved",
    names(comparisons),
    vapply(comparisons, `[[`, "", "id"),
    setup$spec$seed,
    .sg_members(records, attempts)
  )
  names(keys) <- names(comparisons)
  keys
}

# Execution -----------------------------------------------------------------------------
#
# The controller lists the store once, reads every record that it uses once,
# and keeps the payloads. A job that finishes reports the record it used; the
# controller reads that record once. Collection uses these payloads and does
# not list the store again (invariant 3).

.sg_execute <- function(setup, st, backend, retry, progress) {
  # Backends may draw random numbers, for example mirai for its socket
  # names. The caller's stream must not move.
  withr::local_preserve_seed()
  spec <- setup$spec
  spec$retry <- retry
  # One timeout for the controller and every worker.
  spec$claim_timeout <- .claim_timeout()
  datasets <- setup$datasets
  kinds <- c(
    "fits",
    if (!is.null(spec$cap$prepare)) "prepared",
    if (length(spec$cap$comparisons)) "comparisons"
  )
  index <- .store_index(st, kinds)
  fits <- .sg_fit_states(setup, index$fits, retry = retry)
  records <- fits$used
  payloads <- vector("list", nrow(fits))
  known <- which(!is.na(records))
  if (length(known)) {
    files <- .store_file(st, "fits", records[known])
    sizes <- file.size(files)
    payloads[known] <- lapply(seq_along(files), function(k) {
      .read_payload(files[[k]], sizes[[k]])
    })
  }
  attempt_ids <- function(rows) {
    vapply(payloads[rows], function(p) p$attempt$id %||% "", character(1))
  }
  todo <- which(fits$state %in% c("fit", "measure"))
  pending <- tabulate(fits$didx[todo], nrow(datasets))
  has_raw <- length(spec$raw_comparisons) > 0L
  done_comparisons <- index$comparisons

  log <- list(kind = character(), id = character())
  n_log <- 0L
  handle <- NULL
  on.exit(if (!is.null(handle)) handle$close(), add = TRUE)

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
    attempts <- attempt_ids(rows)
    keys <- .sg_compare_keys(setup, datasets$dataset_id[[d]], records[rows], attempts)
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
      attempts = attempts,
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

  # What a finished job changes. A fit reports the record it used, and the
  # controller reads that record here, once.
  complete <- function(job, result) {
    if (job$kind == "prepare") {
      blocked[[job$cidx]] <<- FALSE
      if (length(waiting[[job$cidx]])) {
        queue <<- c(waiting[[job$cidx]], queue[seq.int(qpos, length.out = length(queue) - qpos + 1L)])
        qpos <<- 1L
        waiting[job$cidx] <<- list(NULL)
      }
    } else if (job$kind == "fit") {
      records[[job$i]] <<- result$record
      payloads[[job$i]] <<- .read_payload(.store_file(st, "fits", result$record))
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

  total <- length(front) + length(todo)
  pb <- NULL
  if (progress && total) {
    pb <- cli::cli_progress_bar("Running study", total = total, clear = FALSE)
  }

  inflight <- list()
  deferred <- list()
  submitted <- 0L
  warned <- 0L

  repeat {
    # The backend opens with the first job, so that a run with nothing to do
    # needs no backend.
    job <- if (is.null(handle) || length(inflight) < handle$capacity) next_job()
    if (!is.null(job)) {
      if (is.null(handle)) {
        handle <- backend$open(.make_worker(spec))
      }
      submitted <- submitted + 1L
      id <- as.character(submitted)
      inflight[[id]] <- job
      handle$submit(id, job[!names(job) %in% c("i", "d")])
      next
    }
    if (!length(inflight)) {
      if (!length(deferred)) {
        break
      }
      # Other processes hold the remaining jobs. Use a job's record as soon
      # as it is in the store. Submit a job again when its claim is gone or
      # stale; the worker then finds it done or runs it.
      Sys.sleep(0.05)
      still <- list()
      for (job in deferred) {
        done <- .job_done(job, spec, st)
        if (!is.null(done)) {
          complete(job, done)
        } else if (.claim_busy(st, job$key, spec$claim_timeout)) {
          still[[length(still) + 1L]] <- job
        } else {
          front <- c(front, list(job))
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
      if (identical(value$status, "busy")) {
        deferred[[length(deferred) + 1L]] <- job
        next
      }
      if (isTRUE(value$ran)) {
        n_log <- n_log + 1L
        log$kind[[n_log]] <- job$kind
        log$id[[n_log]] <- job$key
      }
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
  list(
    executed = .sg_frame(log, length(log$kind)),
    records = records,
    payloads = payloads,
    comparisons = done_comparisons
  )
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
#
# Build the results from the records that .sg_execute() used. Comparisons of
# saved measurements are stored under a key of the records they read; a run
# on the same records reads them.

.sg_collect <- function(setup, st, exec, study, saved = TRUE) {
  spec <- setup$spec
  datasets <- setup$datasets
  fits <- setup$fits
  payloads <- exec$payloads
  records <- exec$records
  attempt_ids <- vapply(payloads, function(p) p$attempt$id %||% "", character(1))
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
    for (d in seq_len(nrow(datasets))) {
      rows <- which(fits$didx == d)
      keys <- .sg_compare_keys(
        setup, datasets$dataset_id[[d]], records[rows], attempt_ids[rows]
      )
      for (k in keys) {
        comparisons[[length(comparisons) + 1L]] <- .read_payload(
          .store_file(st, "comparisons", k)
        )$value
      }
    }
  }
  if (saved) {
    plain <- Filter(function(x) !length(x$needs), spec$cap$comparisons)
    keys <- .sg_saved_keys(setup, plain, records, attempt_ids)
    for (name in names(plain)) {
      file <- .store_file(st, "comparisons", keys[[name]])
      value <- if (keys[[name]] %in% exec$comparisons) {
        .read_payload(file)$value
      } else {
        computed <- .sg_compare(
          study$comparisons[name],
          measurements,
          attempts,
          stream = c("compare", spec$seed)
        )[[name]]
        payload <- serialize(list(value = computed), NULL, xdr = FALSE)
        unserialize(.publish_raw(payload, file))$value
      }
      comparisons[[length(comparisons) + 1L]] <- value
    }
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
