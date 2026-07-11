# Controller-only checkpoint scaling benchmark.
#
# Run manually from the package root, for example:
#   Rscript inst/benchmarks/checkpoint-scaling.R 10000 100000
#
# This intentionally does not execute fitters or workers. It measures the
# controller-side checkpoint write/read path using synthetic completed task
# results. Results are printed as a tibble; temporary checkpoints are removed.

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("This benchmark requires devtools.")
}
devtools::load_all(quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
task_counts <- if (length(args)) as.integer(args) else c(10000L, 100000L)

measure_one <- function(n_tasks) {
  checkpoint_root <- tempfile("bayesim-checkpoint-benchmark-")
  init_checkpoint_dir(
    checkpoint_root,
    config_fingerprint = "checkpoint-scaling-benchmark"
  )
  on.exit(unlink(checkpoint_root, recursive = TRUE), add = TRUE)

  ids <- sprintf("d001_f001_r%05d", seq_len(n_tasks))
  task_grid <- tibble::tibble(
    data_idx = 1L,
    fit_idx = 1L,
    rep_idx = seq_len(n_tasks),
    task_id = ids,
    rng_seed = rep(list(1:7), n_tasks),
    status = "success"
  )
  task_results <- lapply(ids, function(id) {
    list(
      task_id = id,
      status = "success",
      metrics = list(stub = list(value = 1))
    )
  })

  gc(reset = TRUE)
  write_time <- system.time(write_checkpoint(
    checkpoint_root,
    task_grid,
    task_results,
    config_fingerprint = "checkpoint-scaling-benchmark"
  ))[["elapsed"]]
  gc_peak <- gc()
  bytes <- sum(
    file.info(list.files(
      checkpoint_root,
      recursive = TRUE,
      full.names = TRUE
    ))$size,
    na.rm = TRUE
  )
  resume_time <- system.time(read_checkpoint(checkpoint_root))[["elapsed"]]

  tibble::tibble(
    tasks = n_tasks,
    controller_wall_seconds = write_time,
    checkpoint_bytes = bytes,
    resume_seconds = resume_time,
    peak_gc_mb = max(gc_peak[, "max used (Mb)"])
  )
}

print(do.call(rbind, lapply(task_counts, measure_one)))
