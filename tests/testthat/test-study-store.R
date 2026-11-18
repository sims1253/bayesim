# The study store: publishing, claims, attempts and what a run reads. The
# invariants these tests check are written at the top of R/study-store.R.

local_store <- function(env = parent.frame()) {
  bayesim:::.store_open(withr::local_tempdir(.local_envir = env))
}

fit_name <- function(i = 1L) strrep(sprintf("%x", i %% 16L), 40L)

# S5 -------------------------------------------------------------------------------

test_that("temporary names hold the host, the process and a nonce, and draw no random numbers", {
  st <- local_store()
  seen <- character()
  writer <- function(x, f) {
    seen <<- c(seen, basename(f))
    writeLines(x, f)
  }
  set.seed(1)
  before <- .Random.seed
  bayesim:::.write_atomic("a", file.path(st$root, "one"), writer)
  bayesim:::.write_atomic("b", file.path(st$root, "two"), writer)
  expect_identical(.Random.seed, before)
  host <- gsub(".", "\\.", bayesim:::.host(), fixed = TRUE)
  pattern <- paste0("^\\.tmp-", host, "-", Sys.getpid(), "-([0-9a-f]{16})-[0-9]+-(one|two)$")
  expect_match(seen, pattern)
  nonce <- sub(pattern, "\\1", seen[[1]])
  # Another process gets another nonce.
  skip_if_not_installed("callr")
  other <- callr::r(function() {
    set.seed(1)
    before <- .Random.seed
    nonce <- bayesim:::.process_nonce()
    list(nonce = nonce, same_seed = identical(.Random.seed, before))
  })
  expect_true(other$same_seed)
  expect_false(identical(other$nonce, nonce))
})

# S3 -------------------------------------------------------------------------------

test_that("opening a store again after an interrupted first open succeeds", {
  root <- withr::local_tempdir()
  # What the first open leaves when it stops before the marker is in place.
  dir.create(file.path(root, "fits"))
  dir.create(file.path(root, "claims"))
  file.create(file.path(root, ".tmp-host-1-0123456789abcdef-1-bayesim-store"))
  st <- bayesim:::.store_open(root)
  expect_true(file.exists(file.path(st$root, "bayesim-store")))
  foreign <- withr::local_tempdir()
  file.create(file.path(foreign, "notes.txt"))
  expect_error(bayesim:::.store_open(foreign), "not empty")
})

# S19 ------------------------------------------------------------------------------

test_that("a store on a file system without hard links is refused when it opens", {
  root <- withr::local_tempdir()
  local_mocked_bindings(.link_file = function(from, to) FALSE)
  expect_error(bayesim:::.store_open(root), "hard links")
  # A write never replaces a file, also when linking fails after the open.
  path <- file.path(root, "one")
  expect_error(bayesim:::.write_atomic("a", path, writeLines), "Could not publish")
  expect_false(file.exists(path))
})

# S6 -------------------------------------------------------------------------------

test_that("a size taken before a record changed does not read as corruption", {
  st <- local_store()
  path <- bayesim:::.store_file(st, "datasets", fit_name())
  bayesim:::.write_payload(list(x = 1), path)
  size <- file.size(path)
  unlink(path)
  bayesim:::.write_payload(list(x = 1:1000), path)
  expect_identical(bayesim:::.read_payload(path, size)$x, 1:1000)
})

# A study of one method whose fit and measurement note each call in `calls`.
noted_study <- function(calls, conditions = "a") {
  note <- function(what) cat(what, "\n", file = calls, append = TRUE, sep = "")
  s <- study("noted", data.frame(condition_id = conditions),
             function(condition, context) rnorm(3))
  s <- with_methods(s, m = study_method(function(data, context) {
    note("fit")
    list(draws = rnorm(5, mean(data)))
  }))
  with_measure(s, "mean", function(artifacts, context) {
    note("mean")
    data.frame(value = mean(artifacts$fit$draws))
  }, needs = "fit")
}

noted <- function(calls, what) {
  if (file.exists(calls)) sum(readLines(calls) == what) else 0L
}

# S29 ------------------------------------------------------------------------------

test_that("two artifact names whose digests share a prefix have two files", {
  st <- local_store()
  # The first eight characters of these names' digests are the same.
  files <- vapply(c("artifact24291", "artifact79363"), function(name) {
    bayesim:::.artifact_file(st, "artifacts", fit_name(), "1", name)
  }, character(1))
  expect_false(files[[1]] == files[[2]])
})

# S1, S2 ---------------------------------------------------------------------------

claimed <- function(x) !is.null(x) && !isFALSE(x)

test_that("a process cannot claim a job again while it holds the claim", {
  st <- local_store()
  first <- bayesim:::.claim(st, fit_name())
  second <- bayesim:::.claim(st, fit_name())
  expect_true(claimed(first))
  expect_false(claimed(second))
})

test_that("of two processes that take over one stale claim, one wins", {
  st <- local_store()
  key <- fit_name()
  dir <- file.path(st$claims, key)
  dir.create(dir)
  old <- "elsewhere_1_0_0"
  file.create(file.path(dir, old))
  # Both contenders find the old claim stale. The other contender takes the
  # claim over after this one has checked it and before this one acts.
  other <- NULL
  entered <- FALSE
  local_mocked_bindings(.claim_stale = function(x, ...) {
    if (!entered) {
      entered <<- TRUE
      other <<- bayesim:::.claim(st, key)
      return(TRUE)
    }
    dir.exists(x) || identical(x, old)
  })
  this <- bayesim:::.claim(st, key)
  expect_equal(claimed(this) + claimed(other), 1)
  expect_length(list.files(dir), 1L)
})

test_that("a process releases only a claim that it owns", {
  st <- local_store()
  key <- fit_name()
  mine <- bayesim:::.claim(st, key)
  # Another process takes the claim over, as it may when this one is slow.
  dir <- file.path(st$claims, key)
  file.rename(file.path(dir, mine), file.path(dir, "elsewhere_1_0_0"))
  bayesim:::.release(st, key, mine)
  expect_identical(list.files(dir), "elsewhere_1_0_0")
})

# A backend that runs each job in its own process, so that no job waits for
# another job that its process has queued. close() lets running jobs run on,
# as mirai daemons without a dispatcher do. `left` keeps those processes, so
# that the test can wait for them and garbage collection does not kill them.
process_backend <- function(left) {
  structure(list(name = "process", open = function(worker) {
    procs <- list()
    list(
      capacity = 2L,
      submit = function(id, job) {
        procs[[id]] <<- callr::r_bg(function(worker, job) worker(job), list(worker, job))
      },
      poll = function() {
        repeat {
          done <- names(procs)[!vapply(procs, function(p) p$is_alive(), logical(1))]
          if (length(done)) break
          Sys.sleep(0.02)
        }
        out <- lapply(done, function(id) {
          list(id = id, value = tryCatch(procs[[id]]$get_result(), error = identity))
        })
        procs[done] <<- NULL
        out
      },
      close = function() {
        for (id in names(procs)) assign(id, procs[[id]], envir = left)
        procs <<- list()
      }
    )
  }), class = "bayesim_backend")
}

test_that("a worker holds its claim until it has published, also after the run stops", {
  skip_on_cran()
  skip_if_not_installed("callr")
  path <- withr::local_tempfile()
  gate <- withr::local_tempdir()
  left <- new.env()
  withr::defer(for (p in as.list(left)) p$kill())
  # A gate that does not open in time fails the test with its name, both in a
  # job and in the test.
  wait_until <- function(ok, what, seconds = 120) {
    deadline <- Sys.time() + seconds
    while (!ok() && Sys.time() < deadline) Sys.sleep(0.02)
    ok()
  }
  wait_for <- function(name) {
    if (!wait_until(function() file.exists(file.path(gate, name)), name)) {
      file.create(file.path(gate, paste0("timeout-", name)))
      stop("gate '", name, "' did not open in time")
    }
  }
  s <- study("hold", data.frame(condition_id = "a"),
             function(condition, context) context$replicate)
  s <- with_methods(s, m = study_method(function(data, context) {
    if (data == 1) {
      file.create(file.path(gate, "started"))
      wait_for("go")
    }
    list(v = data)
  }))
  s <- with_measure(s, "v", function(artifacts, context) {
    if (artifacts$fit$v == 2) {
      wait_for("started")
      stop("bad measure")
    }
    data.frame(value = artifacts$fit$v)
  }, needs = "fit")
  expect_error(run_study(s, 2, path = path, backend = process_backend(left)), "bad measure")
  claims_dir <- file.path(path, "claims")
  claims <- function() list.files(claims_dir, pattern = "^[^.]")
  records <- function() list.files(file.path(path, "fits"), recursive = TRUE, full.names = TRUE)
  # The run has stopped. The fit of replicate 1 still runs: it holds its
  # claim and has published nothing.
  expect_length(claims(), 1L)
  expect_length(records(), 0L)
  file.create(file.path(gate, "go"))
  if (!wait_until(function() !length(claims()), "release")) {
    fail("the worker did not release its claim within 120 s of 'go'")
  }
  expect_length(list.files(gate, pattern = "^timeout-"), 0L)
  record <- records()
  expect_length(record, 1L)
  if (length(record) == 1L) {
    expect_true(bayesim:::.is_success_record(basename(sub("[.]rds$", "", record))))
    # Linking the record to its name changes its status time. Removing the
    # claim directory is the last change to claims/. The record came first.
    expect_lte(
      as.numeric(file.info(record)$ctime),
      as.numeric(file.info(claims_dir)$mtime)
    )
  }
})

# Invariant 1 ------------------------------------------------------------------------
#
# Two runs of one study start at once. The claim timeout is zero, so each run
# takes over the other's claims and both do every job. Each callback waits
# until the other run has reached it too, so both generate every dataset and
# fit every method before either publishes. The data and the fits depend on
# the process ID, so the two attempts at each job differ.

race_script <- function(path, gate, peers, extra = FALSE) {
  library(bayesim)
  options(bayesim.claim_timeout = 0, race.peers = peers)
  meet <- function(step) {
    file.create(file.path(gate, paste0(step, "-", Sys.getpid())))
    deadline <- Sys.time() + 30
    while (length(list.files(gate, pattern = paste0("^", step, "-"))) <
           getOption("race.peers") && Sys.time() < deadline) {
      Sys.sleep(0.01)
    }
  }
  s <- study("race", data.frame(condition_id = "a"), function(condition, context) {
    meet(paste0("generate-", context$dataset_id))
    Sys.getpid() + context$replicate / 10
  })
  fit <- function(data, context) {
    meet(paste0("fit-", context$dataset_id, "-", context$method_id))
    list(data = data, own = Sys.getpid() + stats::runif(1))
  }
  s <- with_methods(s, a = study_method(fit), b = study_method(fit))
  s <- with_measure(s, "own", function(artifacts, context) {
    data.frame(value = artifacts$fit$own)
  }, needs = "fit")
  s <- with_measure(s, "data", function(artifacts, context) {
    data.frame(value = artifacts$fit$data)
  }, needs = "fit")
  s <- with_comparison(s, "gap", function(results, artifacts, context) {
    data.frame(value = artifacts$a$fit$own - artifacts$b$fit$own)
  }, needs = "fit")
  s <- with_retention(s, c("data", "fit"))
  if (extra) {
    # Measurements on what the store kept: the dataset and the fit.
    s <- with_measure(s, "stored", function(artifacts, context) {
      data.frame(value = artifacts$data)
    }, needs = "data")
    s <- with_measure(s, "kept", function(artifacts, context) {
      data.frame(value = artifacts$fit$own)
    }, needs = "fit")
  }
  run_study(s, 3, path = path, progress = FALSE)
}

test_that("two runs of a non-deterministic study at once leave a consistent store", {
  skip_on_cran()
  skip_if_not_installed("callr")
  path <- withr::local_tempfile()
  gate <- withr::local_tempdir()
  ps <- lapply(1:2, function(i) callr::r_bg(race_script, list(path, gate, 2L)))
  for (p in ps) p$wait(300000)
  runs <- lapply(ps, function(p) p$get_result())
  # Both runs fitted every method on every dataset.
  expect_length(list.files(gate, pattern = "^fit-"), 2L * 6L)
  # Each run reports the records that it used, and both used the store's.
  expect_identical(runs[[1]]$measurements, runs[[2]]$measurements)
  expect_identical(runs[[1]]$comparisons, runs[[2]]$comparisons)
  m <- runs[[1]]$measurements
  value <- function(d, method, what) {
    m$value[m$dataset_id == d & m$method_id == method & m$measurement == what]
  }
  gap <- runs[[1]]$comparisons
  for (d in unique(m$dataset_id)) {
    # The comparison saw the fits whose measurements the run reports.
    expect_equal(gap$value[gap$dataset_id == d], value(d, "a", "own") - value(d, "b", "own"))
    # Both methods saw one dataset.
    expect_identical(value(d, "a", "data"), value(d, "b", "data"))
  }
  # What the store kept is what the measurements saw. This run measures the
  # kept data and fits and fits nothing.
  check <- callr::r(race_script, list(path, gate, 1L, TRUE))
  expect_length(list.files(gate, pattern = "^fit-"), 2L * 6L)
  k <- check$measurements
  for (what in c("own", "data")) {
    expect_identical(k[k$measurement == what, "value"], m[m$measurement == what, "value"])
  }
  expect_identical(k[k$measurement == "kept", "value"], k[k$measurement == "own", "value"])
  expect_identical(k[k$measurement == "stored", "value"], k[k$measurement == "data", "value"])
})

# S4 ---------------------------------------------------------------------------------

test_that("a failed fit that another process retries during a run does not break the run", {
  skip_on_cran()
  skip_if_not_installed("callr")
  path <- withr::local_tempfile()
  flag <- withr::local_tempfile()
  writeLines("fail", flag)
  base <- study("retry", data.frame(condition_id = "a"), function(condition, context) 1)
  base <- with_methods(
    base,
    flaky = study_method(function(data, context) {
      if (readLines(flag) == "fail") stop("no convergence")
      list(v = 1)
    }),
    steady = study_method(function(data, context) list(v = 2))
  )
  base <- with_measure(base, "v", function(artifacts, context) {
    data.frame(value = artifacts$fit$v)
  }, needs = "fit")
  base <- with_comparison(base, "count", function(results, artifacts, context) {
    data.frame(value = length(artifacts))
  }, needs = "fit")
  base <- with_retention(base, "fit")
  first <- run_study(base, 1, path = path)
  expect_equal(first$attempts$status, c("failed", "success"))
  # A new method's fit retries the failed fit in another process, after this
  # run has planned and before it runs the comparison.
  writeLines("ok", flag)
  s <- with_methods(base, trigger = study_method(function(data, context) {
    callr::r(function(s, path) {
      bayesim::run_study(s, 1, path = path, retry = TRUE)
      NULL
    }, list(base, path))
    list(v = 3)
  }))
  run <- run_study(s, 1, path = path)
  expect_equal(run$attempts$status, c("failed", "success", "success"))
  expect_equal(run$comparisons$value, 2)
  # The next run uses the fit that succeeded.
  later <- run_study(s, 1, path = path)
  expect_equal(later$attempts$status, rep("success", 3))
  expect_equal(later$comparisons$value, 3)
})

# S10 --------------------------------------------------------------------------------

test_that("an artifact that was discarded and is needed again is computed again", {
  path <- withr::local_tempfile()
  calls <- withr::local_tempfile()
  s <- noted_study(calls)
  first <- run_study(s, 2, path = path)
  more <- with_measure(s, "max", function(artifacts, context) {
    data.frame(value = max(artifacts$fit$draws))
  }, needs = "fit")
  again <- run_study(more, 2, path = path)
  # The fits run again; the old measurement is not computed again.
  expect_equal(noted(calls, "fit"), 4)
  expect_equal(noted(calls, "mean"), 2)
  expect_identical(again$attempts, first$attempts)
  cmp <- with_comparison(more, "spread", function(results, artifacts, context) {
    data.frame(value = stats::sd(artifacts$m$fit$draws))
  }, needs = "fit")
  compared <- run_study(cmp, 2, path = path)
  expect_equal(noted(calls, "fit"), 6)
  # The seeds make these the same fits as on a new store.
  fresh <- run_study(cmp, 2)
  expect_identical(again$measurements, fresh$measurements)
  expect_identical(compared$comparisons, fresh$comparisons)
})

# S25 --------------------------------------------------------------------------------

test_that("a run that newly keeps a fit keeps the transient copy before its record names it", {
  path <- withr::local_tempfile()
  flag <- withr::local_tempfile()
  s <- with_comparison(noted_study(withr::local_tempfile()), "spread", function(results, artifacts, context) {
    if (file.exists(flag)) stop("bad comparison")
    data.frame(value = stats::sd(artifacts$m$fit$draws))
  }, needs = "fit")
  # The comparison fails, so the fits stay in the transient files.
  file.create(flag)
  expect_error(run_study(s, 1, path = path), "bad comparison")
  unlink(flag)
  expect_length(list.files(file.path(path, "transient"), recursive = TRUE), 1L)
  run_study(with_retention(s, "fit"), 1, path = path)
  # The record that keeps the fit names a kept file.
  expect_length(list.files(file.path(path, "artifacts"), recursive = TRUE), 1L)
  expect_length(list.files(file.path(path, "transient"), recursive = TRUE), 0L)
})

# S26 --------------------------------------------------------------------------------

test_that("a process that loses the race to publish a transient artifact uses its own value", {
  st <- local_store()
  spec <- list(retention = character(), compress = FALSE)
  # Another process published the file first, and a comparison deleted it.
  local_mocked_bindings(.write_raw = function(...) FALSE)
  out <- bayesim:::.job_publish_artifacts(st, spec, fit_name(), "1", list(fit = 42))
  expect_identical(out$artifacts, list(fit = 42))
  expect_length(out$files, 0L)
})

# S30 --------------------------------------------------------------------------------

# A backend that runs each job in this process, inside `wrap`.
inline_backend <- function(wrap = function(run, job) run()) {
  structure(list(name = "inline", open = function(worker) {
    done <- list()
    list(
      capacity = 1L,
      submit = function(id, job) {
        done[[id]] <<- list(id = id, value = wrap(function() worker(job), job))
      },
      poll = function() {
        out <- unname(done)
        done <<- list()
        out
      },
      close = function() done <<- list()
    )
  }), class = "bayesim_backend")
}

# Stop a run that waits on a claim more than `limit` times.
local_wait_limit <- function(limit, first = NULL, env = parent.frame()) {
  busy <- bayesim:::.claim_busy
  n <- 0L
  local_mocked_bindings(.claim_busy = function(...) {
    n <<- n + 1L
    if (n > limit) stop("the run kept waiting on a claim")
    if (n == 1L && !is.null(first)) first()
    busy(...)
  }, .package = "bayesim", .env = env)
}

test_that("a run uses a busy job's record as soon as another process has published it", {
  path <- withr::local_tempfile()
  other <- withr::local_tempfile()
  s <- noted_study(withr::local_tempfile())
  run_study(s, 1, path = other)
  # Another process holds every claim and keeps it. It publishes its records
  # while this run waits.
  local_mocked_bindings(.claim = function(...) NULL, .claim_busy = function(...) TRUE)
  local_wait_limit(100L, first = function() {
    dirs <- list.dirs(other, recursive = FALSE)
    file.copy(dirs[basename(dirs) != "claims"], path, recursive = TRUE)
  })
  run <- run_study(s, 1, path = path, backend = inline_backend())
  expect_equal(nrow(run$executed), 0L)
  expect_identical(run$measurements, run_study(s, 1, path = other)$measurements)
})

test_that("a worker judges a claim with the timeout of the process that runs the study", {
  path <- withr::local_tempfile()
  s <- noted_study(withr::local_tempfile())
  # A claim of another host that is ten seconds old is on every job.
  old <- sprintf("elsewhere_1_%.0f_x", (as.numeric(Sys.time()) - 10) * 1000)
  claim_old <- function(run, job) {
    dir <- file.path(path, "claims", job$key)
    if (!dir.exists(dir)) {
      dir.create(dir)
      file.create(file.path(dir, old))
    }
    # The worker's own option says that the claim is live.
    withr::with_options(list(bayesim.claim_timeout = 3600), run())
  }
  withr::local_options(bayesim.claim_timeout = 1)
  local_wait_limit(50L)
  run <- run_study(s, 1, path = path, backend = inline_backend(claim_old))
  # The study has one job, and the worker ran it.
  expect_identical(run$executed$kind, "fit")
})

# S11 --------------------------------------------------------------------------------

# Method `first` changes the dataset in place; method `second` reads it.
mutating_study <- function(make, change) {
  s <- study("shared", data.frame(condition_id = "a"),
             function(condition, context) make())
  s <- with_methods(
    s,
    first = study_method(function(data, context) {
      change(data)
      list(x = data$x)
    }),
    second = study_method(function(data, context) list(x = data$x))
  )
  with_measure(s, "x", function(artifacts, context) {
    data.frame(value = artifacts$fit$x)
  }, needs = "fit")
}

seen_by <- function(run, method) {
  run$measurements$value[run$measurements$method_id == method]
}

test_that("a job does not see what another job did to an environment dataset", {
  s <- mutating_study(
    function() {
      e <- new.env()
      e$x <- 1
      e
    },
    function(data) data$x <- 99
  )
  for (kept in list(character(), "data")) {
    run <- run_study(with_retention(s, kept), 1)
    expect_equal(seen_by(run, "first"), 99)
    expect_equal(seen_by(run, "second"), 1)
  }
})

test_that("a job does not see what another job did to a data.table dataset", {
  skip_if_not_installed("data.table")
  s <- mutating_study(
    function() data.table::data.table(x = 1),
    function(data) data.table::set(data, j = "x", value = 99)
  )
  for (kept in list(character(), "data")) {
    run <- run_study(with_retention(s, kept), 1)
    expect_equal(seen_by(run, "first"), 99)
    expect_equal(seen_by(run, "second"), 1)
  }
})

test_that("a job does not see what another job did to the prepared value", {
  s <- study(
    "prepared",
    data.frame(condition_id = "a"),
    function(condition, context) 1,
    prepare = function(condition, context) {
      e <- new.env()
      e$x <- 1
      e
    }
  )
  s <- with_methods(
    s,
    first = study_method(function(data, context) {
      context$prepared$x <- 99
      list(x = context$prepared$x)
    }),
    second = study_method(function(data, context) list(x = context$prepared$x))
  )
  s <- with_measure(s, "x", function(artifacts, context) {
    data.frame(value = artifacts$fit$x)
  }, needs = "fit")
  run <- run_study(s, 1)
  expect_equal(seen_by(run, "second"), 1)
})

test_that("a job does not see what another job did to a captured environment", {
  state <- new.env()
  state$n <- 0
  s <- study("captured", data.frame(condition_id = "a"), function(condition, context) 1)
  # The fit reads `state` and changes it through a local name. An assignment
  # to `state$n` itself would make codetools treat `state` as a local
  # variable, and capture would not copy it.
  s <- with_methods(s, m = study_method(function(data, context) {
    counter <- state
    counter$n <- counter$n + 1
    list(n = counter$n)
  }))
  s <- with_measure(s, "n", function(artifacts, context) {
    data.frame(value = artifacts$fit$n)
  }, needs = "fit")
  run <- run_study(s, 3)
  expect_equal(run$measurements$value, c(1, 1, 1))
})

# S12 --------------------------------------------------------------------------------

test_that("a run that newly keeps the data writes the datasets before its records", {
  path <- withr::local_tempfile()
  flag <- withr::local_tempfile()
  s <- study("keep", data.frame(condition_id = "a"), function(condition, context) {
    if (file.exists(flag)) stop("no data")
    rnorm(2)
  })
  s <- with_methods(s, m = study_method(function(data, context) list(v = data[[1]])))
  s <- with_measure(s, "v", function(artifacts, context) {
    data.frame(value = artifacts$fit$v)
  }, needs = "fit")
  run_study(s, 2, path = path)
  kept <- run_study(with_retention(s, "data"), 2, path = path)
  datasets <- list.files(file.path(path, "datasets"), recursive = TRUE)
  expect_length(datasets, 2L)
  # When the data cannot be made, no record that promises them is written.
  path <- withr::local_tempfile()
  run_study(s, 1, path = path)
  file.create(flag)
  expect_error(run_study(with_retention(s, "data"), 1, path = path), "no data")
  expect_length(list.files(file.path(path, "fits"), recursive = TRUE), 1L)
})

# S13 --------------------------------------------------------------------------------

test_that("a run reuses every saved measurement of a fit, whichever record holds it", {
  path <- withr::local_tempfile()
  calls <- withr::local_tempfile()
  s <- with_retention(noted_study(calls), "fit")
  s$measures <- list()
  note <- function(what) cat(what, "\n", file = calls, append = TRUE, sep = "")
  one <- function(artifacts, context) {
    note("one")
    data.frame(value = 1)
  }
  two <- function(artifacts, context) {
    note("two")
    data.frame(value = 2)
  }
  run_study(with_measure(s, "one", one, needs = "fit"), 1, path = path)
  run_study(with_measure(s, "two", two, needs = "fit"), 1, path = path)
  unlink(calls)
  both <- with_measure(with_measure(s, "one", one, needs = "fit"), "two", two, needs = "fit")
  run <- run_study(both, 1, path = path)
  expect_equal(nrow(run$measurements), 2)
  expect_equal(noted(calls, "fit") + noted(calls, "one") + noted(calls, "two"), 0)
})

# S14 --------------------------------------------------------------------------------

test_that("a run with nothing to do opens no backend", {
  skip_on_cran()
  skip_if_not_installed("mirai")
  path <- withr::local_tempfile()
  s <- noted_study(withr::local_tempfile())
  run_study(s, 1, path = path)
  run <- run_study(s, 1, path = path, backend = mirai_backend())
  expect_equal(nrow(run$executed), 0)
})

test_that("comparisons of saved measurements are stored under the records they read", {
  path <- withr::local_tempfile()
  calls <- withr::local_tempfile()
  s <- with_comparison(noted_study(calls), "count", function(results, artifacts, context) {
    cat("count\n", file = calls, append = TRUE)
    data.frame(value = nrow(results))
  }, by = character())
  first <- run_study(s, 2, path = path)
  expect_equal(noted(calls, "count"), 1)
  second <- run_study(s, 2, path = path)
  expect_equal(noted(calls, "count"), 1)
  expect_identical(second$comparisons, first$comparisons)
  third <- run_study(s, 3, path = path)
  expect_equal(noted(calls, "count"), 2)
  expect_equal(third$comparisons$value, 3)
})
