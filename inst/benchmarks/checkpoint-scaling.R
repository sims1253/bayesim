# Controller-only checkpoint scaling benchmark.
#
# Run manually from the package root, for example:
#   Rscript inst/benchmarks/checkpoint-scaling.R 1000 10000
# Optional second-style arguments are task counts; add "--batch N" via env
# BAYESIM_BENCH_BATCH to change the checkpoint cadence (default 100).
#
# This intentionally does not execute fitters or workers. It measures the
# controller-side checkpoint write/read path exactly the way execute_tasks()
# drives it: canonical outcomes built with new_task_result() are pushed
# through a filesystem RunStore in replicate batches, so each checkpoint
# commit appends only the new outcomes as an immutable shard. Outcomes
# constructed as plain lists are silently dropped by write_checkpoint()
# (it filters for bayesim_task_result objects), which is why this benchmark
# must go through the constructor. Results are printed as a tibble;
# temporary checkpoints are removed.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("This benchmark requires devtools.")
}
devtools::load_all(quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
task_counts <- if (length(args)) as.integer(args) else c(1000L, 10000L)
batch_size <- as.integer(Sys.getenv("BAYESIM_BENCH_BATCH", "100"))

FINGERPRINT <- "checkpoint-scaling-benchmark"

# Small but realistic metric payload: the flattened fields a typical
# posterior_summary + coverage study produces per task.
stub_metrics <- function(i) {
  list(
    posterior_summary__mean__x = 0.5 + i * 1e-6,
    posterior_summary__sd__x = 0.1,
    posterior_summary__mean__sigma = 1.0,
    coverage__by_param__x = 1
  )
}

measure_one <- function(n_tasks) {
  checkpoint_root <- tempfile("bayesim-checkpoint-benchmark-")
  init_checkpoint_dir(
    checkpoint_root,
    config_fingerprint = FINGERPRINT
  )
  on.exit(unlink(checkpoint_root, recursive = TRUE), add = TRUE)

  ids <- sprintf("d001_f001_r%05d", seq_len(n_tasks))
  task_grid <- tibble::tibble(
    data_idx = 1L,
    fit_idx = 1L,
    rep_idx = seq_len(n_tasks),
    task_id = ids,
    rng_seed = rep(list(1:7), n_tasks),
    status = "pending"
  )
  outcomes <- lapply(seq_len(n_tasks), function(i) {
    new_task_result(
      task_id = ids[[i]],
      status = "success",
      metrics = stub_metrics(i),
      timing = list(total = 0.25)
    )
  })
  names(outcomes) <- ids

  # The same store the controller uses: each write() appends one immutable
  # outcome shard containing only the outcomes new since the previous commit.
  store <- new_run_store(
    result_path = checkpoint_root,
    config_fingerprint = FINGERPRINT,
    keep_checkpoints = 2L
  )

  gc(reset = TRUE)
  write_time <- system.time({
    for (batch_start in seq.int(1L, n_tasks, by = batch_size)) {
      batch_end <- min(batch_start + batch_size - 1L, n_tasks)
      batch_idx <- batch_start:batch_end
      task_grid$status[batch_idx] <- "success"
      # Mirror execute_tasks(): pass every non-NULL outcome accumulated so
      # far; the delta store writes only the ones that changed.
      store$write(
        task_grid = task_grid,
        task_results = outcomes[seq_len(batch_end)]
      )
    }
  })[["elapsed"]]
  gc_peak <- gc()

  outcome_files <- list.files(
    file.path(checkpoint_root, "outcomes"),
    pattern = "\\.rds$"
  )
  bytes <- sum(
    file.info(list.files(
      checkpoint_root,
      recursive = TRUE,
      full.names = TRUE
    ))$size,
    na.rm = TRUE
  )
  resume_time <- system.time(read_checkpoint(checkpoint_root))[["elapsed"]]

  # Integrity check: every outcome must have landed in the store.
  checkpoint <- read_checkpoint(checkpoint_root)
  n_persisted <- length(checkpoint$task_outcomes)

  # gc() column naming differs across R versions ("max used (Mb)" pre-4.6 vs
  # paired "max used"/"(Mb)" columns); the last column is the Mb value of
  # "max used" in both layouts.
  peak_gc_mb <- if ("max used (Mb)" %in% colnames(gc_peak)) {
    max(gc_peak[, "max used (Mb)"])
  } else {
    max(gc_peak[, ncol(gc_peak)])
  }

  tibble::tibble(
    tasks = n_tasks,
    batch_size = batch_size,
    controller_wall_seconds = write_time,
    outcome_shards = length(outcome_files),
    outcomes_persisted = n_persisted,
    checkpoint_bytes = bytes,
    resume_seconds = resume_time,
    peak_gc_mb = peak_gc_mb
  )
}

results <- do.call(rbind, lapply(task_counts, measure_one))
if (any(results$outcomes_persisted != results$tasks)) {
  warning(
    "Benchmark integrity failure: persisted outcome count does not match ",
    "task count. The checkpoint path is dropping outcomes."
  )
}
print(results)
