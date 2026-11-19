# Backends -------------------------------------------------------------------------
#
# A backend decides where jobs run. run_study() needs only four operations
# from it, so a new backend (a cluster scheduler, for example) is a small
# function.

#' Study backends
#'
#' A backend runs the jobs of [run_study()]. `sequential_backend()` runs them
#' one at a time in this R process. `mirai_backend()` runs them in parallel
#' on [mirai](https://mirai.r-lib.org) daemons.
#'
#' With `mirai_backend()`, the study keeps at most two jobs per worker
#' submitted at a time, and submits the next job as soon as one finishes.
#' There is no barrier between batches, so a slow fit holds up one worker
#' only. Each worker writes its record to the store when its job finishes.
#'
#' @section Backend contract:
#' A backend is a list of class `bayesim_backend` with a `name` and a
#' function `open(worker)`. `worker` is a function of one job; it returns a
#' short list and writes all results to the store itself. `open()` returns a
#' list of:
#'
#' * `capacity`: the number of jobs that may be submitted at a time.
#' * `submit(id, job)`: start `worker(job)`. `id` is a string.
#' * `poll()`: wait until at least one submitted job finishes, then return a
#'   list of `list(id = , value = )` for the finished jobs. `value` is the
#'   value of `worker(job)`, or an error object when the job could not run.
#' * `close()`: stop unfinished jobs and free resources. [run_study()] calls
#'   it once, also after an error or an interrupt.
#'
#' A backend must send `worker` to each process once, not with each job: it
#' contains all the study's functions.
#'
#' @param n Number of daemons to start. With `NULL`, use the daemons that are
#'   already set with [mirai::daemons()].
#' @param ... Further arguments for [mirai::daemons()] when `n` is given.
#'
#' @return A `bayesim_backend`.
#' @name study-backends
NULL

#' @rdname study-backends
#' @export
sequential_backend <- function() {
  structure(
    list(
      name = "sequential",
      open = function(worker) {
        queue <- list()
        list(
          capacity = 1L,
          submit = function(id, job) {
            queue[[id]] <<- job
          },
          poll = function() {
            id <- names(queue)[[1L]]
            job <- queue[[id]]
            queue[[id]] <<- NULL
            list(list(id = id, value = worker(job)))
          },
          close = function() {
            queue <<- list()
          }
        )
      }
    ),
    class = "bayesim_backend"
  )
}

#' @rdname study-backends
#' @export
mirai_backend <- function(n = NULL, ...) {
  if (!is.null(n)) {
    .sg_count(n, "n")
  }
  dots <- list(...)
  structure(
    list(
      name = "mirai",
      n = n,
      open = function(worker) {
        if (is.null(n)) {
          profile <- "default"
          info <- mirai::status(.compute = profile)
          workers <- info$connections
          if (!isTRUE(workers > 0L)) {
            cli::cli_abort(c(
              "No mirai daemons are running.",
              i = "Call {.code mirai::daemons(n)} first, or use {.code mirai_backend(n)}."
            ))
          }
          own <- FALSE
        } else {
          profile <- paste0("bayesim-", Sys.getpid(), "-", .next_tmp_id())
          do.call(
            mirai::daemons,
            c(list(n = n, dispatcher = TRUE, .compute = profile), dots)
          )
          workers <- n
          own <- TRUE
        }
        ok <- FALSE
        on.exit(if (!ok && own) mirai::daemons(0L, .compute = profile))
        # Send the worker once. Later jobs refer to it by name.
        mirai::everywhere(
          {
            NULL
          },
          .bayesim_worker = worker,
          .compute = profile
        )
        # The worker exists only on the daemons, so the call is quoted.
        call <- quote(.bayesim_worker(job))
        tasks <- list()
        ok <- TRUE
        list(
          capacity = 2L * workers,
          submit = function(id, job) {
            tasks[[id]] <<- mirai::mirai(
              .expr = call,
              job = job,
              .compute = profile
            )
          },
          poll = function() {
            mirai::race_mirai(tasks, .compute = profile)
            done <- names(tasks)[!vapply(tasks, mirai::unresolved, logical(1))]
            out <- lapply(done, function(id) {
              list(id = id, value = tasks[[id]]$data)
            })
            tasks[done] <<- NULL
            out
          },
          close = function() {
            for (task in tasks) mirai::stop_mirai(task)
            tasks <<- list()
            if (own) {
              mirai::daemons(0L, .compute = profile)
            }
          }
        )
      }
    ),
    class = "bayesim_backend"
  )
}

#' @export
print.bayesim_backend <- function(x, ...) {
  cat("<bayesim backend> ", x$name,
      if (!is.null(x$n)) paste0(" (", x$n, " workers)"), "\n", sep = "")
  invisible(x)
}
