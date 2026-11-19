# Study engine benchmark.
#
# Run manually from the package root, on a machine that is otherwise idle:
#   Rscript inst/benchmarks/study-engine.R
#
# The acceptance tests in tests/testthat/test-study-acceptance.R check what
# the engine does: how many records a resume reads, how much store work each
# job costs, and that a slow fit holds up only its own worker. They assert no
# time limit, so that they give the same verdict on any machine. This script
# measures how fast the engine does that work, against four limits:
#
#   resume    a run with nothing left to do, on a store of 5,000 records,
#             takes under 1 s (the fastest of 3 runs)
#   overhead  the engine's CPU time per job, beyond the user's functions, is
#             under 2 ms at 500 and at 5,000 jobs, and the 5,000-job figure is
#             under 1.5 times the 500-job figure plus 0.25 ms
#   growth    the total CPU time per job at 5,000 jobs is under 3 times the
#             figure at 500 jobs
#   straggle  40 fits on 4 workers, three of them ten times slower than the
#             rest, take under 1.2 times the ideal schedule
#
# The script prints one row per measurement and exits with status 1 if any
# limit is exceeded. On a loaded machine the figures say little: compare runs
# on the same idle machine before and after a change.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("This benchmark requires devtools.")
}
devtools::load_all(quiet = TRUE)

timing_study <- function() {
  s <- study("timing", data.frame(condition_id = "a", n = 5L),
             function(condition, context) rnorm(condition$n))
  s <- with_methods(s, m = study_method(function(data, context) list(mu = mean(data))))
  with_measure(s, "mu", function(artifacts, context) data.frame(value = artifacts$fit$mu),
               needs = "fit")
}

cpu <- function(t) t[["user.self"]] + t[["sys.self"]]

bench_resume <- function(n = 5000) {
  path <- tempfile()
  on.exit(unlink(path, recursive = TRUE))
  s <- timing_study()
  run_study(s, n, path = path, progress = FALSE)
  elapsed <- vapply(1:3, function(i) {
    system.time(run_study(s, n, path = path, progress = FALSE))[["elapsed"]]
  }, numeric(1))
  value <- min(elapsed)
  data.frame(
    measurement = "resume",
    value = sprintf("%.3f s for %d records", value, n),
    limit = "under 1 s",
    pass = value < 1
  )
}

# Overhead is the CPU time of run_study() minus the CPU time that the same
# user functions take in a plain loop, divided by the number of jobs.
bench_overhead <- function(sizes = c(500, 5000)) {
  s <- timing_study()
  user_cpu <- function(n) {
    cpu(system.time(for (i in seq_len(n)) {
      data <- s$generate(s$conditions[1, ], list())
      fit <- s$methods$m$fit(data, list())
      s$measures$mu$compute(list(fit = fit), list())
    }))
  }
  ms <- vapply(sizes, function(n) {
    min(vapply(1:2, function(i) {
      path <- tempfile()
      on.exit(unlink(path, recursive = TRUE))
      engine <- cpu(system.time(run_study(s, n, path = path, progress = FALSE)))
      1000 * (engine - user_cpu(n)) / n
    }, numeric(1)))
  }, numeric(1))
  data.frame(
    measurement = "overhead",
    value = sprintf("%.2f ms per job at %d, %.2f ms at %d", ms[1], sizes[1], ms[2], sizes[2]),
    limit = "under 2 ms; second under 1.5 x first + 0.25",
    pass = max(ms) < 2 && ms[2] < 1.5 * ms[1] + 0.25
  )
}

# Store operations do not show work in memory that grows faster than the
# number of jobs, for example a vector copied on each job. CPU time does, so
# this measurement compares the CPU time per job at 500 and 5,000 jobs. The
# sequential backend runs everything in this process, and CPU time leaves out
# the time that other processes take, so machine load moves the two sizes
# together. On a 12-core workstation at a load average of about 24, three
# measurements of the ratio gave 0.98, 1.07 and 1.09. Work that grows with the
# square of the number of jobs makes the ratio approach 10 at these sizes. The
# limit of 3 lies between the two. Each size runs twice, and the faster run
# counts. The acceptance tests used to assert this limit; they now assert
# counts only, which give the same verdict on any machine.
bench_growth <- function(sizes = c(500, 5000)) {
  s <- timing_study()
  per_job <- vapply(sizes, function(n) {
    min(vapply(1:2, function(i) {
      path <- tempfile()
      on.exit(unlink(path, recursive = TRUE))
      cpu(system.time(run_study(s, n, path = path, progress = FALSE))) / n
    }, numeric(1)))
  }, numeric(1))
  ratio <- per_job[2] / per_job[1]
  data.frame(
    measurement = "growth",
    value = sprintf("CPU time per job at %d is %.2f x the time at %d", sizes[2], ratio, sizes[1]),
    limit = "under 3 x",
    pass = ratio < 3
  )
}

bench_straggle <- function(unit = 0.2, workers = 4) {
  duration <- rep(unit, 40)
  duration[c(5, 15, 25)] <- 10 * unit
  s <- study("straggle", data.frame(condition_id = "a"),
             function(condition, context) context$replicate)
  s <- with_methods(s, m = study_method(function(data, context) {
    Sys.sleep(duration[[data]])
    list()
  }))
  mirai::daemons(workers, .compute = "default")
  on.exit(mirai::daemons(0, .compute = "default"))
  # Warm the daemons, so that only scheduling counts.
  run_study(s, 1, backend = mirai_backend(), progress = FALSE)
  wall <- system.time(
    run_study(s, length(duration), backend = mirai_backend(), progress = FALSE)
  )[["elapsed"]]
  ideal <- max(sum(duration) / workers, max(duration))
  data.frame(
    measurement = "straggle",
    value = sprintf("%.2f s against an ideal of %.2f s (%.2fx)", wall, ideal, wall / ideal),
    limit = "under 1.2 x ideal",
    pass = wall < 1.2 * ideal
  )
}

results <- rbind(bench_resume(), bench_overhead(), bench_growth(), bench_straggle())
print(results, right = FALSE, row.names = FALSE)
if (!all(results$pass)) {
  quit(status = 1)
}
