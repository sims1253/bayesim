# Acceptance tests for the study engine. None of them asserts a time limit;
# inst/benchmarks/study-engine.R measures the engine's speed.

# A study whose measurement reports the data and the method's own draw, so
# that tests can see both random streams.
acceptance_study <- function(methods = c("a", "b"), conditions = c("x", "y")) {
  s <- study(
    "acceptance",
    data.frame(condition_id = conditions, n = 4L),
    function(condition, context) rnorm(condition$n)
  )
  m <- study_method(function(data, context) list(data = data, own = runif(1)))
  s <- do.call(with_methods, c(list(s), stats::setNames(rep(list(m), length(methods)), methods)))
  with_measure(
    s,
    "trace",
    function(artifacts, context) {
      data.frame(
        first = artifacts$fit$data[[1L]],
        own = artifacts$fit$own,
        extra = runif(1)
      )
    },
    needs = "fit"
  )
}

sorted_rows <- function(x) {
  x <- x[order(x$dataset_id, x$method_id), , drop = FALSE]
  rownames(x) <- NULL
  x
}

# 1. Pairing ------------------------------------------------------------------

test_that("data depend on seed, condition and replicate; method draws also on the method", {
  two <- run_study(acceptance_study(c("a", "b")), 3, seed = 7)$measurements
  one <- run_study(acceptance_study("b", c("y", "x")), 3, seed = 7)$measurements
  # Same datasets, whatever the methods and the condition order.
  b2 <- sorted_rows(two[two$method_id == "b", ])
  b1 <- sorted_rows(one)
  expect_identical(b1$dataset_id, b2$dataset_id)
  expect_identical(b1$first, b2$first)
  # Both methods see the same data, but draw their own numbers.
  a2 <- sorted_rows(two[two$method_id == "a", ])
  expect_identical(a2$first, b2$first)
  expect_false(any(a2$own == b2$own))
  # A method's stream does not depend on the other methods.
  expect_identical(b1$own, b2$own)
  expect_identical(b1$extra, b2$extra)
  # Another seed gives other data.
  other <- run_study(acceptance_study("b"), 3, seed = 8)$measurements
  expect_false(any(other$first %in% b2$first))
})

# 2. Invariance -----------------------------------------------------------------

test_that("results do not depend on the backend, the workers or the store", {
  skip_on_cran()
  s <- acceptance_study()
  ref <- run_study(s, 6)
  for (n in c(1L, 3L)) {
    run <- run_study(s, 6, backend = mirai_backend(n))
    expect_identical(run$measurements, ref$measurements)
    expect_identical(run$attempts[setdiff(names(run$attempts), "elapsed")],
                     ref$attempts[setdiff(names(ref$attempts), "elapsed")])
  }
  # Built in two steps in one store: the same rows.
  path <- withr::local_tempfile()
  run_study(s, 2, path = path)
  expect_identical(run_study(s, 6, path = path)$measurements, ref$measurements)
})

# 3. Extension -------------------------------------------------------------------

test_that("new replicates, methods and measurements run only the new work", {
  path <- withr::local_tempfile()
  s <- acceptance_study("a")
  first <- run_study(s, 2, path = path)
  expect_equal(nrow(first$executed), 4)

  more <- run_study(s, 3, path = path)
  expect_equal(nrow(more$executed), 2)
  expect_identical(more$measurements[more$measurements$replicate <= 2, ], first$measurements)

  wider <- acceptance_study(c("a", "b"))
  w <- run_study(wider, 3, path = path)
  expect_equal(nrow(w$executed), 6)
  expect_identical(
    sorted_rows(w$measurements[w$measurements$method_id == "a", ]),
    sorted_rows(more$measurements)
  )

  # A new measurement on a kept fit measures again; it does not fit again.
  measured <- with_measure(
    with_retention(wider, "fit"),
    "size",
    function(artifacts, context) data.frame(value = length(artifacts$fit$data)),
    needs = "fit"
  )
  path2 <- withr::local_tempfile()
  base <- run_study(with_retention(wider, "fit"), 2, path = path2)
  again <- run_study(measured, 2, path = path2)
  expect_equal(nrow(again$executed), 8)
  expect_identical(
    sorted_rows(again$measurements[again$measurements$measurement == "trace", ])$own,
    sorted_rows(base$measurements)$own
  )
  expect_identical(again$attempts$elapsed, base$attempts$elapsed)
  expect_equal(nrow(run_study(measured, 2, path = path2)$executed), 0)
})

# 4. Fingerprint -------------------------------------------------------------------

test_that("a change to anything a function uses changes the key", {
  path <- withr::local_tempfile()
  make <- function(shift, helper, version = "1") {
    # Capture refuses a promise that is not forced.
    force(shift)
    force(helper)
    s <- study("fingerprint", data.frame(condition_id = "a"),
               function(condition, context) rnorm(3))
    s <- with_methods(s, m = study_method(
      function(data, context) list(v = helper(mean(data)) + shift),
      version = version
    ))
    with_measure(s, "v", function(artifacts, context) data.frame(value = artifacts$fit$v),
                 needs = "fit")
  }
  add <- function(x) x + 1
  base <- run_study(make(0, add), 2, path = path)
  expect_equal(nrow(base$executed), 2)
  expect_equal(nrow(run_study(make(0, add), 2, path = path)$executed), 0)
  expect_equal(nrow(run_study(make(1, add), 2, path = path)$executed), 2)
  times <- function(x) x * 2
  expect_equal(nrow(run_study(make(0, times), 2, path = path)$executed), 2)
  expect_equal(nrow(run_study(make(0, add, "2"), 2, path = path)$executed), 2)
  # Every variant has its own records; none overwrote another.
  expect_equal(nrow(run_study(make(0, add), 2, path = path)$executed), 0)

  # Values that cannot reach a worker fail at plan time, by name.
  con <- file(withr::local_tempfile(), "w")
  on.exit(close(con), add = TRUE)
  bad <- study("bad", data.frame(condition_id = "a"), function(condition, context) 1)
  bad <- with_methods(bad, slow = study_method(function(data, context) {
    writeLines("x", con)
    list()
  }))
  expect_error(plan_study(bad, 1), "slow.*captures.*con.*cannot be sent")
  big <- runif(2e5)
  huge <- with_methods(
    study("big", data.frame(condition_id = "a"), function(condition, context) 1),
    heavy = study_method(function(data, context) list(m = mean(big)))
  )
  withr::local_options(bayesim.capture_limit = 1e5)
  expect_error(plan_study(huge, 1), "heavy.*captures.*big.*bytes")
})

# 5. No-op resume --------------------------------------------------------------------
#
# Tests 5 to 7 assert what the engine does, not how long it takes, so that
# they give the same verdict on any machine. The timings live in
# inst/benchmarks/study-engine.R.

timing_study <- function(calls = NULL) {
  # Capture refuses a promise that is not forced.
  force(calls)
  note <- function(stage) {
    if (!is.null(calls)) cat(stage, "\n", file = calls, append = TRUE, sep = "")
  }
  s <- study("timing", data.frame(condition_id = "a", n = 5L),
             function(condition, context) {
               note("generate")
               rnorm(condition$n)
             })
  s <- with_methods(s, m = study_method(function(data, context) {
    note("fit")
    list(mu = mean(data))
  }))
  with_measure(s, "mu", function(artifacts, context) {
    note("measure")
    data.frame(value = artifacts$fit$mu)
  }, needs = "fit")
}

store_ops <- function(expr) bayesim:::.with_store_trace(expr)

test_that("a run with nothing left to do walks the fits once, reads each record once and runs nothing", {
  n <- 500
  path <- withr::local_tempfile()
  calls <- withr::local_tempfile()
  s <- timing_study(calls)
  run_study(s, n, path = path)
  expect_equal(as.vector(table(readLines(calls))[c("generate", "fit", "measure")]), rep(n, 3))
  unlink(calls)

  traced <- store_ops(run_study(s, n, path = path))
  run <- traced$value
  ops <- traced$ops
  expect_equal(nrow(run$executed), 0)
  expect_equal(nrow(run$measurements), n)
  # No user function ran: none of them wrote to the log.
  expect_false(file.exists(calls))
  # The only store work is one walk of the fits, the only kind of record
  # this study needs, and one read per record. Nothing is claimed, probed,
  # written or deleted, and collection lists nothing again.
  expect_setequal(unique(ops$op), c("walk", "read"))
  expect_equal(sum(ops$op == "walk"), 1)
  reads <- ops$path[ops$op == "read"]
  records <- list.files(file.path(path, "fits"), recursive = TRUE, full.names = TRUE)
  expect_length(records, n)
  expect_equal(anyDuplicated(reads), 0L)
  expect_setequal(normalizePath(reads), normalizePath(records))
})

# 6. Work per job ----------------------------------------------------------------------
#
# Store operations are counted exactly. For n jobs on a new store, the
# worker of each job claims it once and releases it once. It probes five
# times. Before it claims the job and again after, it probes for the job's
# record and for failed attempts (one listing of the record's directory);
# then it probes for a stored dataset. It writes its record, and
# the controller reads that record once. The run as a whole writes the store
# marker and walks the fits once, when it plans; it deletes nothing. Any
# operation that grows faster than n breaks these equalities at one of the
# sizes.

test_that("each job claims, releases, probes five times, writes and reads once", {
  ops_of <- c("claim", "release", "probe", "write", "read", "walk", "delete")
  counts <- function(n) {
    path <- withr::local_tempfile()
    ops <- store_ops(run_study(timing_study(), n, path = path))$ops
    expect_true(all(ops$op %in% ops_of))
    vapply(ops_of, function(op) sum(ops$op == op), numeric(1))
  }
  for (n in c(50, 500)) {
    expect_equal(
      counts(n),
      c(claim = n, release = n, probe = 5 * n, write = n + 1, read = n, walk = 1, delete = 0),
      label = paste(n, "jobs")
    )
  }
})

# 7. Stragglers ------------------------------------------------------------------------
#
# The engine submits fits in study order, so the three slow fits (replicates
# 1 to 3) go first and take three of the four workers. Each slow fit waits
# until all 37 fast fits have left their mark, then records how many it saw.
# If the engine kept the fast fits back until the slow ones finished, a slow
# fit would see fewer than 37. The deadline of two minutes only stops a
# broken engine from hanging the suite; a working one releases the slow fits
# as soon as the last fast fit finishes.

test_that("slow fits hold up one worker, not the run", {
  skip_on_cran()
  marks <- withr::local_tempdir()
  n <- 40
  slow <- 1:3
  s <- study("straggle", data.frame(condition_id = "a"),
             function(condition, context) context$replicate)
  s <- with_methods(s, m = study_method(function(data, context) {
    if (data %in% slow) {
      deadline <- Sys.time() + 120
      seen <- function() length(list.files(marks, pattern = "^fast-"))
      while (seen() < n - length(slow) && Sys.time() < deadline) Sys.sleep(0.02)
      writeLines(as.character(seen()), file.path(marks, paste0("slow-", data)))
    } else {
      file.create(file.path(marks, paste0("fast-", data)))
    }
    list()
  }))
  mirai::daemons(4, .compute = "default")
  on.exit(mirai::daemons(0, .compute = "default"), add = TRUE)
  run <- run_study(s, n, backend = mirai_backend())
  expect_equal(nrow(run$executed), n)
  expect_length(list.files(marks, pattern = "^fast-"), n - length(slow))
  seen <- vapply(paste0("slow-", slow), function(f) {
    as.integer(readLines(file.path(marks, f)))
  }, integer(1))
  expect_equal(unname(seen), rep(n - length(slow), length(slow)))
})

# 8. Failure ----------------------------------------------------------------------------

test_that("a failed fit is recorded, kept, and fitted again only on request", {
  path <- withr::local_tempfile()
  flag <- withr::local_tempfile()
  writeLines("fail", flag)
  s <- study("failure", data.frame(condition_id = "a"), function(condition, context) 1)
  s <- with_methods(
    s,
    flaky = study_method(function(data, context) {
      if (readLines(flag) == "fail") stop("no convergence")
      list(v = 1)
    }),
    steady = study_method(function(data, context) list(v = 2))
  )
  s <- with_measure(s, "v", function(artifacts, context) data.frame(value = artifacts$fit$v),
                    needs = "fit")
  first <- run_study(s, 2, path = path)
  failed <- first$attempts[first$attempts$method_id == "flaky", ]
  expect_equal(failed$status, c("failed", "failed"))
  expect_match(failed$error, "no convergence")
  expect_setequal(first$measurements$method_id, "steady")

  writeLines("ok", flag)
  expect_equal(nrow(run_study(s, 2, path = path)$executed), 0)
  again <- run_study(s, 2, path = path, retry = TRUE)
  expect_equal(nrow(again$executed), 2)
  expect_equal(again$attempts$status, rep("success", 4))

  # An error outside `fit` stops the run, says where, and keeps earlier records.
  path <- withr::local_tempfile()
  broken <- with_measure(s, "broken", function(artifacts, context) {
    if (context$replicate == 2) stop("bad measure")
    data.frame(value = 0)
  }, needs = "fit")
  err <- expect_error(run_study(broken, 3, path = path), class = "bayesim_study_error")
  expect_match(conditionMessage(err), "bad measure")
  expect_match(cli::ansi_strip(conditionMessage(err)), "Stage measure, condition a, replicate 2")
  st <- bayesim:::.store_open(path)
  expect_equal(length(bayesim:::.store_index(st, "fits")$fits), 2)
})

# 9. Interrupt ------------------------------------------------------------------------

test_that("an interrupted run leaves a consistent store and resumes", {
  skip_on_cran()
  skip_if_not_installed("callr")
  path <- withr::local_tempfile()
  started <- withr::local_tempfile()
  script <- function(path, started, backend) {
    library(bayesim)
    s <- study("interrupt", data.frame(condition_id = "a"), function(condition, context) 1)
    s <- with_methods(s, m = study_method(function(data, context) {
      cat("x\n", file = started, append = TRUE)
      Sys.sleep(0.2)
      list(v = 1)
    }))
    s <- with_measure(s, "v", function(artifacts, context) data.frame(value = 1), needs = "fit")
    b <- if (backend == "mirai") mirai_backend(2) else sequential_backend()
    run_study(s, 200, path = path, backend = b)
  }
  for (backend in c("sequential", "mirai")) {
    unlink(c(path, started), recursive = TRUE)
    p <- callr::r_bg(script, list(path, started, backend))
    deadline <- Sys.time() + 60
    while (Sys.time() < deadline &&
           (!file.exists(started) || length(readLines(started)) < 6)) {
      p$wait(100)
    }
    p$interrupt()
    p$wait(30000)
    expect_false(p$is_alive())
    st <- bayesim:::.store_open(path)
    done <- bayesim:::.store_index(st, "fits")$fits
    expect_gt(length(done), 0)
    expect_lt(length(done), 200)
    # No partial files, and no claim held by a live process.
    expect_length(
      list.files(path, pattern = "^\\.tmp-", recursive = TRUE, all.files = TRUE, include.dirs = TRUE),
      0
    )
    tokens <- list.files(st$claims, recursive = TRUE)
    expect_true(all(vapply(basename(tokens), bayesim:::.claim_stale, logical(1))))
    # The rest runs; finished fits are not run again.
    rest <- callr::r(script, list(path, started, "sequential"))
    expect_equal(nrow(rest$executed), 200 - length(done))
    expect_equal(nrow(rest$measurements), 200)
  }
})

# 10. Concurrent writers -----------------------------------------------------------------

test_that("several processes share one store and run each fit once", {
  skip_on_cran()
  skip_if_not_installed("callr")
  path <- withr::local_tempfile()
  log <- withr::local_tempfile()
  script <- function(path, log) {
    library(bayesim)
    s <- study("concurrent", data.frame(condition_id = c("a", "b")),
               function(condition, context) rnorm(2))
    s <- with_methods(s, m = study_method(function(data, context) {
      cat(context$dataset_id, "\n", file = log, append = TRUE)
      Sys.sleep(0.01)
      list(v = data[[1L]])
    }), n = study_method(function(data, context) list(v = runif(1))))
    s <- with_measure(s, "v", function(artifacts, context) data.frame(value = artifacts$fit$v),
                      needs = "fit")
    run_study(s, 30, path = path)
  }
  ps <- lapply(1:3, function(i) callr::r_bg(script, list(path, log)))
  for (p in ps) p$wait(120000)
  runs <- lapply(ps, function(p) p$get_result())
  lines <- trimws(readLines(log))
  expect_length(lines, 60)
  expect_false(anyDuplicated(lines) > 0)
  expect_equal(sum(vapply(runs, function(r) nrow(r$executed), 0)), 120)
  for (r in runs[-1]) expect_identical(r$measurements, runs[[1]]$measurements)
  fresh <- callr::r(script, list(withr::local_tempfile(), withr::local_tempfile()))
  expect_identical(fresh$measurements, runs[[1]]$measurements)
})

# 11. Calibration -------------------------------------------------------------------------

test_that("the conjugate model is calibrated and a wrong prior is not", {
  skip_on_cran()
  s <- study(
    "calibration",
    data.frame(condition_id = "n20", n = 20L),
    prior_draws(
      prior = function(condition) c(b1 = rnorm(1), b2 = rnorm(1)),
      simulate = function(truth, condition) {
        X <- cbind(b1 = 1, b2 = rnorm(condition$n))
        list(X = X, y = drop(X %*% truth) + rnorm(condition$n))
      }
    )
  )
  s <- with_methods(
    s,
    exact = conjugate_lm(draws = 200),
    narrow = conjugate_lm(prior_sd = 0.1, draws = 200)
  )
  s <- with_measure(s, "coverage", measure_coverage(0.9))
  s <- with_measure(s, "bias", measure_bias())
  s <- with_measure(s, "rmse", measure_rmse())
  s <- with_measure(s, "rank", measure_sbc_rank())
  s <- with_comparison(s, "rmse_gain", compare_paired("rmse", reference = "exact"))
  run <- run_study(s, 400)
  sm <- summarise_study(run)
  exact <- sm[sm$method_id == "exact", ]
  narrow <- sm[sm$method_id == "narrow", ]

  cov <- exact[exact$measurement == "coverage", ]
  expect_true(all(abs(cov$estimate - 0.9) < 4 * cov$mcse))
  rank <- exact[exact$measurement == "rank", ]
  expect_true(all(abs(rank$estimate - 0.5) < 4 * rank$mcse))
  b <- exact[exact$measurement == "bias", ]
  expect_true(all(abs(b$estimate) < 4 * b$mcse))
  # The narrow prior shrinks to zero; its intervals miss the truth.
  ncov <- narrow[narrow$measurement == "coverage", ]
  expect_true(all(ncov$estimate < 0.9 - 4 * ncov$mcse))

  cmp <- summarise_study(run, "comparisons", target_mcse = 0.01)
  expect_setequal(cmp$target, c("b1", "b2"))
  expect_true(all(cmp$estimate > 4 * cmp$mcse))
  expect_true(all(c("n_needed", "done") %in% names(cmp)))
  # Pairing helps: the paired MCSE is below that of two separate means.
  m <- run$measurements[run$measurements$measurement == "rmse", ]
  for (t in c("b1", "b2")) {
    a <- m$value[m$method_id == "exact" & m$target == t]
    z <- m$value[m$method_id == "narrow" & m$target == t]
    unpaired <- sqrt(var(a) / length(a) + var(z) / length(z))
    expect_lt(cmp$mcse[cmp$target == t], unpaired)
  }
})

# Vignettes --------------------------------------------------------------------------------

# The likelihood and primed-priors vignettes do not run their brms fits
# (`eval = FALSE`), so they render without brms or CmdStan too.
test_that("the grammar vignettes render against the package", {
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("callr")
  skip_if_not(rmarkdown::pandoc_available())
  root <- test_path("..", "..", "vignettes")
  skip_if_not(dir.exists(root), "vignette sources are not installed")
  names <- c(
    "experimental-study-grammar",
    "experimental-method-comparison",
    "experimental-likelihood-study",
    "experimental-primed-priors"
  )
  for (name in names) {
    out <- withr::local_tempdir()
    file.copy(file.path(root, paste0(name, ".Rmd")), out)
    html <- callr::r(
      function(input, out) rmarkdown::render(input, output_dir = out, quiet = TRUE,
                                             envir = new.env()),
      list(file.path(out, paste0(name, ".Rmd")), out)
    )
    expect_true(file.exists(html))
  }
})
