test_that("existing result paths are never silently reused", {
  gen <- function(data_spec, task_ctx) {
    list(
      train = data.frame(y = rnorm(8), x = rnorm(8)),
      test = NULL,
      response = "y",
      true_params = c(beta = 1),
      vars_of_interest = "beta",
      meta = list()
    )
  }
  path <- file.path(withr::local_tempdir(), "run")
  cfg <- simulation_config(
    data_grid = data.frame(n = 1L),
    fit_grid = data.frame(model = "lm"),
    data_generator = gen,
    fitter = LinearRegressionFitter(n_draws = 10L),
    metrics = list(),
    n_replicates = 1L,
    seed = 1L,
    result_path = path
  )
  run_simulation(cfg, resume = "never", progress = FALSE)

  incompatible <- cfg
  incompatible@data_grid <- data.frame(n = 2L)
  expect_error(
    run_simulation(incompatible, resume = "auto", progress = FALSE),
    class = "bayesim_checkpoint_error"
  )
})

test_that("lifecycle predicates partition every declared task status", {
  expect_identical(
    is_terminal_task_status(TASK_STATUSES),
    TASK_STATUSES %in% TASK_TERMINAL_STATUSES
  )
  expect_identical(
    is_resumable_task_status(TASK_STATUSES),
    TASK_STATUSES %in% TASK_RESUMABLE_STATUSES
  )
  expect_true(all(
    is_terminal_task_status(TASK_STATUSES) !=
      is_resumable_task_status(TASK_STATUSES)
  ))
})

test_that("zero error budget stops at the first error and resume runs placeholders", {
  gen <- function(data_spec, task_ctx) {
    if (task_ctx$rep_idx == 2L) {
      stop("intentional integration failure")
    }
    list(
      train = data.frame(y = rnorm(8), x = rnorm(8)),
      test = NULL,
      response = "y",
      true_params = c(beta = 1),
      vars_of_interest = "beta",
      meta = list()
    )
  }
  path <- file.path(withr::local_tempdir(), "run")
  common <- list(
    data_grid = data.frame(n = 1L),
    fit_grid = data.frame(model = "lm"),
    data_generator = gen,
    fitter = LinearRegressionFitter(n_draws = 10L),
    metrics = list(),
    n_replicates = 3L,
    seed = 1L,
    result_path = path,
    checkpoint_every = 1L
  )
  first <- do.call(
    simulation_config,
    c(common, list(max_errors = 0))
  )
  first_result <- run_simulation(first, resume = "never", progress = FALSE)
  expect_equal(sum(first_result$summary$status == "skipped"), 1L)

  continued <- do.call(
    simulation_config,
    c(common, list(max_errors = Inf))
  )
  continued_result <- run_simulation(
    continued,
    resume = "auto",
    progress = FALSE
  )
  expect_equal(sum(continued_result$summary$status == "skipped"), 0L)
  expect_equal(sum(continued_result$summary$status == "failed"), 1L)
})

test_that("uninterrupted, resumed, and two-worker runs have lifecycle parity", {
  gen <- function(data_spec, task_ctx) {
    if (task_ctx$rep_idx == 2L) {
      stop("intentional parity failure")
    }
    list(
      train = data.frame(y = rnorm(12), x = rnorm(12)),
      test = NULL,
      response = "y",
      true_params = c(beta = 1),
      vars_of_interest = "beta",
      meta = list()
    )
  }
  make_config <- function(path, max_errors = Inf) {
    simulation_config(
      data_grid = data.frame(n = 12L),
      fit_grid = data.frame(model = "lm"),
      data_generator = gen,
      fitter = LinearRegressionFitter(n_draws = 10L),
      metrics = list(),
      n_replicates = 4L,
      seed = 1204L,
      result_path = path,
      checkpoint_every = 1L,
      max_errors = max_errors
    )
  }
  run_quietly <- function(config, resume, workers) {
    run_simulation(
      config,
      resume = resume,
      workers = workers,
      progress = FALSE,
      verbose = FALSE
    )
  }
  canonical_summary <- function(result) {
    timing_columns <- grep("time|timing", names(result$summary), value = TRUE)
    result$summary[
      order(result$summary$task_id),
      setdiff(
        names(result$summary),
        timing_columns
      ),
      drop = FALSE
    ]
  }

  uninterrupted <- run_quietly(
    make_config(file.path(withr::local_tempdir(), "uninterrupted")),
    resume = "never",
    workers = 1L
  )

  resume_path <- file.path(withr::local_tempdir(), "resumed")
  interrupted <- run_quietly(
    make_config(resume_path, max_errors = 0),
    resume = "never",
    workers = 1L
  )
  expect_equal(sum(interrupted$summary$status == "skipped"), 2L)
  resumed <- run_quietly(
    make_config(resume_path),
    resume = "auto",
    workers = 1L
  )

  parallel <- run_quietly(
    make_config(file.path(withr::local_tempdir(), "parallel")),
    resume = "never",
    workers = 2L
  )

  expect_equal(canonical_summary(resumed), canonical_summary(uninterrupted))
  expect_equal(canonical_summary(parallel), canonical_summary(uninterrupted))
})

test_that("checkpoint round-trip retains canonical task truth", {
  gen <- function(data_spec, task_ctx) {
    list(
      train = data.frame(y = rnorm(8), x = rnorm(8)),
      test = NULL,
      response = "y",
      true_params = c(beta = 1.25),
      vars_of_interest = "beta",
      meta = list()
    )
  }
  path <- file.path(withr::local_tempdir(), "run")
  cfg <- simulation_config(
    data_grid = data.frame(n = 1L),
    fit_grid = data.frame(model = "lm"),
    data_generator = gen,
    fitter = LinearRegressionFitter(n_draws = 10L),
    metrics = list(),
    n_replicates = 1L,
    seed = 1L,
    result_path = path
  )
  run_simulation(cfg, resume = "never", progress = FALSE)
  checkpoint <- read_checkpoint(path)

  expect_equal(checkpoint$task_outcomes[[1]]$truth[["beta"]], 1.25)
  expect_false("truth" %in% names(checkpoint$task_outcomes[[1]]$metrics))
})

test_that("memory and filesystem run stores round-trip across several batches", {
  # Three writes of five outcomes each: the second and third writes append a
  # batch to an existing store, matching how execute_tasks() checkpoints
  # successive batches. The in-memory adapter must accumulate (not re-flatten
  # from scratch) and read() must return the same outcomes and the same flat
  # results_df as the filesystem adapter.
  make_outcome <- function(i) {
    new_task_result(
      task_id = sprintf("d001_f001_r%05d", i),
      status = "success",
      metrics = list(bias = 0.25 + i),
      diagnostics = list(rhat = 1.01),
      timing = list(total = 0.125 + i),
      warnings = "test warning",
      truth = c(beta = 1.25)
    )
  }
  batches <- lapply(
    list(1:5, 6:10, 11:15),
    function(ids) lapply(ids, make_outcome)
  )
  # The engine always checkpoints the FULL task grid with updated statuses
  # (never a growing grid), so each write passes all 15 rows.
  make_grid <- function(n_success) {
    data.frame(
      task_id = sprintf("d001_f001_r%05d", seq_len(15L)),
      status = ifelse(seq_len(15L) <= n_success, "success", "pending"),
      stop_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  path <- file.path(withr::local_tempdir(), "run")

  memory <- new_run_store()
  memory$initialize()
  filesystem <- new_run_store(
    result_path = path,
    config_fingerprint = "store-test",
    checkpoint_format = "rds"
  )
  filesystem$initialize()

  completed <- c(5L, 10L, 15L)
  for (b in seq_along(batches)) {
    memory$write(make_grid(completed[[b]]), batches[[b]])
    filesystem$write(make_grid(completed[[b]]), batches[[b]])
  }

  memory_checkpoint <- memory$read()
  filesystem_checkpoint <- filesystem$read()

  expect_length(memory_checkpoint$task_outcomes, 15L)
  expect_length(filesystem_checkpoint$task_outcomes, 15L)

  # The adapters agree on the accumulated outcomes and the derived flat view.
  expect_equal(
    lapply(memory_checkpoint$task_outcomes, `[[`, "task_id"),
    lapply(filesystem_checkpoint$task_outcomes, `[[`, "task_id")
  )
  for (field in c("status", "metrics", "diagnostics", "warnings", "truth")) {
    expect_equal(
      lapply(memory_checkpoint$task_outcomes, `[[`, field),
      lapply(filesystem_checkpoint$task_outcomes, `[[`, field)
    )
  }
  expect_equal(memory_checkpoint$results_df, filesystem_checkpoint$results_df)
  # The flat view matches a direct flattening of the accumulated outcomes.
  expect_equal(
    memory_checkpoint$results_df,
    results_to_dataframe(filesystem_checkpoint$task_outcomes)
  )
})

test_that("resume rejects retention widening for completed outcomes", {
  gen <- function(data_spec, task_ctx) {
    list(
      train = data.frame(y = rnorm(12), x = rnorm(12)),
      test = NULL,
      response = "y",
      true_params = c(beta = 1),
      vars_of_interest = "beta",
      meta = list()
    )
  }
  path <- file.path(withr::local_tempdir(), "retention-run")
  cfg <- simulation_config(
    data_grid = data.frame(n = 12L),
    fit_grid = data.frame(model = "lm"),
    data_generator = gen,
    fitter = LinearRegressionFitter(n_draws = 10L),
    metrics = list(),
    n_replicates = 2L,
    seed = 5L,
    result_path = path,
    retain = "minimal"
  )
  run_simulation(cfg, resume = "never", progress = FALSE, verbose = FALSE)

  widened <- cfg
  widened@retain <- resolve_retention_spec(c("metrics", "draws"))
  expect_error(
    run_simulation(
      widened,
      resume = "auto",
      progress = FALSE,
      verbose = FALSE
    ),
    regexp = "retention.*draws|draws.*retention",
    class = "bayesim_checkpoint_error"
  )
})

test_that("resume error budget includes persisted failed outcomes", {
  gen <- function(data_spec, task_ctx) stop("persistent failure")
  make_config <- function(path, max_errors) {
    simulation_config(
      data_grid = data.frame(n = 1L),
      fit_grid = data.frame(model = "failure"),
      data_generator = gen,
      fitter = LinearRegressionFitter(n_draws = 10L),
      metrics = list(),
      n_replicates = 3L,
      seed = 8L,
      result_path = path,
      checkpoint_every = 1L,
      max_errors = max_errors
    )
  }

  unchanged_path <- file.path(withr::local_tempdir(), "unchanged-budget")
  original <- make_config(unchanged_path, 1)
  first <- run_simulation(
    original,
    resume = "never",
    progress = FALSE,
    verbose = FALSE
  )
  expect_equal(
    as.integer(table(first$summary$status)[c("failed", "skipped")]),
    c(1, 2)
  )

  unchanged <- run_simulation(
    make_config(unchanged_path, 1),
    resume = "auto",
    progress = FALSE,
    verbose = FALSE
  )
  expect_equal(
    as.integer(table(unchanged$summary$status)[c("failed", "skipped")]),
    c(1, 2)
  )
  # The carried-over failure count exhausts the budget before any task runs,
  # so the restored placeholders must be re-labeled as policy stops instead of
  # staying pending with no recorded reason (#64).
  expect_equal(
    unchanged$task_grid$status,
    c("failed", "skipped", "skipped")
  )
  expect_equal(
    unchanged$task_grid$stop_reason,
    c(NA_character_, "max_errors", "max_errors")
  )

  raised_path <- file.path(withr::local_tempdir(), "raised-budget")
  run_simulation(
    make_config(raised_path, 1),
    resume = "never",
    progress = FALSE,
    verbose = FALSE
  )
  raised <- run_simulation(
    make_config(raised_path, 2),
    resume = "auto",
    progress = FALSE,
    verbose = FALSE
  )
  expect_equal(
    as.integer(table(raised$summary$status)[c("failed", "skipped")]),
    c(2, 1)
  )
})

test_that("filesystem checkpoints append one-outcome shards", {
  gen <- function(data_spec, task_ctx) {
    list(
      train = data.frame(y = rnorm(10), x = rnorm(10)),
      test = NULL,
      response = "y",
      true_params = c(beta = 1),
      vars_of_interest = "beta",
      meta = list()
    )
  }
  path <- file.path(withr::local_tempdir(), "sharded-run")
  cfg <- simulation_config(
    data_grid = data.frame(n = 10L),
    fit_grid = data.frame(model = "lm"),
    data_generator = gen,
    fitter = LinearRegressionFitter(n_draws = 10L),
    metrics = list(),
    n_replicates = 4L,
    seed = 9L,
    result_path = path,
    checkpoint_every = 1L,
    keep_checkpoints = 2L
  )
  run_simulation(cfg, resume = "never", progress = FALSE, verbose = FALSE)

  shard_paths <- list.files(
    file.path(path, "outcomes"),
    pattern = "^shard_[0-9]{6}\\.rds$",
    full.names = TRUE
  )
  expect_equal(length(shard_paths), 4L)
  expect_true(all(
    vapply(shard_paths, function(x) length(readRDS(x)), integer(1)) == 1L
  ))

  checkpoint <- read_checkpoint(path)
  expect_equal(length(checkpoint$meta$outcome_shards %||% list()), 0L)
  expect_equal(checkpoint$meta$storage_mode, "delta-v1")
  expect_equal(length(checkpoint$task_outcomes), 4L)

  # A single corrupt historical shard is recovered from its immutable mirror,
  # even after checkpoint pruning has removed the snapshot that created it.
  historical_shard <- shard_paths[[2L]]
  writeBin(charToRaw("corrupt"), historical_shard)
  expect_warning(
    recovered <- read_checkpoint(path),
    "Recovered a corrupt shard"
  )
  expect_equal(length(recovered$task_outcomes), 4L)
  expect_equal(length(list_checkpoints(path)), 2L)

  # Descriptor corruption also fails over to an independent mirror instead of
  # silently omitting a committed outcome shard.
  historical_mirror <- sub("\\.rds$", ".mirror.rds", historical_shard)
  expect_true(file.copy(historical_mirror, historical_shard, overwrite = TRUE))
  historical_descriptor <- sub("\\.rds$", ".json", historical_shard)
  writeBin(charToRaw("corrupt"), historical_descriptor)
  expect_warning(
    recovered <- read_checkpoint(path),
    "Recovered a corrupt shard descriptor"
  )
  expect_equal(length(recovered$task_outcomes), 4L)
})

test_that("live filesystem RunStore writes checksum only each new shard", {
  path <- file.path(withr::local_tempdir(), "linear-store")
  store <- new_run_store(
    result_path = path,
    config_fingerprint = "linear-store-test",
    checkpoint_format = "rds",
    keep_checkpoints = 2L
  )
  store$initialize()

  shard_checks <- 0L
  original_checksum <- compute_file_checksum
  local_mocked_bindings(
    compute_file_checksum = function(path) {
      if (
        identical(basename(dirname(path)), "outcomes") &&
          !grepl("\\.mirror\\.rds$", path)
      ) {
        shard_checks <<- shard_checks + 1L
      }
      original_checksum(path)
    },
    .package = "bayesim"
  )

  outcomes <- list()
  for (i in seq_len(8L)) {
    outcome <- new_task_result(
      task_id = sprintf("task-%02d", i),
      status = "success",
      metrics = list(value = i),
      timing = list(total = 0)
    )
    outcomes <- c(outcomes, list(outcome))
    grid <- data.frame(
      task_id = vapply(outcomes, function(x) x$task_id, character(1)),
      status = "success",
      stop_reason = NA_character_
    )
    store$write(grid, outcomes)
  }

  # One checksum for each newly appended shard. Rechecking accumulated shards
  # would produce 8 + 7 + ... + 1 checks instead.
  expect_equal(shard_checks, 8L)
})

test_that("delta RunStore serializes O(total tasks) ledger rows", {
  path <- file.path(withr::local_tempdir(), "delta-scaling")
  store <- new_run_store(
    result_path = path,
    config_fingerprint = "delta-scaling-test",
    checkpoint_format = "rds",
    keep_checkpoints = 2L
  )
  store$initialize()

  ledger_rows <- 0L
  original_write <- write_redundant_shard
  local_mocked_bindings(
    write_redundant_shard = function(object, path) {
      if (identical(basename(dirname(path)), "ledger")) {
        ledger_rows <<- ledger_rows + nrow(object)
      }
      original_write(object, path)
    },
    .package = "bayesim"
  )

  n <- 20L
  outcomes <- list()
  grid <- data.frame(
    task_id = sprintf("task-%02d", seq_len(n)),
    status = "pending",
    stop_reason = NA_character_
  )
  for (i in seq_len(n)) {
    outcomes[[i]] <- new_task_result(
      task_id = grid$task_id[[i]],
      status = "success",
      metrics = list(value = i),
      timing = list(total = 0)
    )
    grid$status[[i]] <- "success"
    store$write(grid, outcomes)
  }

  # One N-row immutable base plus one changed row per later checkpoint.
  expect_lte(ledger_rows, 2L * n)
  latest_meta <- jsonlite::read_json(file.path(
    path,
    "checkpoints",
    "cp_000020",
    "meta.json"
  ))
  # Use exact extraction: `$outcome_shards` partially matches the singular
  # `outcome_shard` field in R lists.
  expect_length(latest_meta[["outcome_shards"]], 0L)
  expect_equal(latest_meta$outcome_shard$n_outcomes, 1L)
})

test_that("configless resume uses the latest effective run policy", {
  failure_path <- file.path(withr::local_tempdir(), "policy-errors")
  failing <- simulation_config(
    data_grid = data.frame(n = 0L, beta = 1, sigma = 1),
    fit_grid = data.frame(model = "invalid-empty-data"),
    data_generator = bayesim_example_data_generator,
    fitter = LinearRegressionFitter(n_draws = 10L),
    metrics = list(),
    n_replicates = 3L,
    seed = 14L,
    result_path = failure_path,
    checkpoint_every = 1L,
    keep_checkpoints = 2L,
    retain = "debug",
    max_errors = 1,
    stop_on = list(
      estimand = "sigma",
      measure = "bias",
      target_mcse = 0.25,
      min_reps = 2L,
      check_every = 2L
    )
  )
  first <- run_simulation(
    failing,
    resume = "never",
    progress = FALSE,
    verbose = FALSE
  )
  expect_equal(
    as.integer(table(first$summary$status)[c("failed", "skipped")]),
    c(1, 2)
  )

  changed <- failing
  changed@max_errors <- 2
  changed@retain <- resolve_retention_spec("minimal")
  second <- resume_simulation(
    failure_path,
    config = changed,
    progress = FALSE,
    verbose = FALSE
  )
  expect_equal(
    as.integer(table(second$summary$status)[c("failed", "skipped")]),
    c(2, 1)
  )

  resumed_failure <- resume_simulation(
    failure_path,
    progress = FALSE,
    verbose = FALSE
  )
  expect_equal(
    as.integer(table(resumed_failure$summary$status)[c("failed", "skipped")]),
    c(2, 1)
  )
  rehydrated <- rehydrate_config_from_manifest(failure_path)
  expect_equal(rehydrated@max_errors, 2)
  expect_equal(rehydrated@stop_on, failing@stop_on)
  expect_equal(rehydrated@checkpoint_every, 1L)
  expect_equal(rehydrated@keep_checkpoints, 2L)
  expect_equal(rehydrated@retain, changed@retain)
  widened <- changed
  widened@retain <- resolve_retention_spec("debug")
  expect_error(
    resume_simulation(
      failure_path,
      config = widened,
      progress = FALSE,
      verbose = FALSE
    ),
    "Cannot widen retention"
  )
})

test_that("a fatal mid-batch error checkpoints the batch's successful outcomes", {
  gen <- function(data_spec, task_ctx) {
    list(
      train = data.frame(y = rnorm(10), x = rnorm(10)),
      test = NULL,
      response = "y",
      true_params = c(beta = 1),
      vars_of_interest = "beta",
      meta = list()
    )
  }
  # Fails fatally (bayesim-classified, non-recoverable) on replicate 2 only;
  # its batch siblings succeed.
  FatalMidBatchFitter <- S7::new_class(
    "FatalMidBatchFitter",
    parent = LinearRegressionFitter,
    properties = list(
      name = S7::new_property(S7::class_character, default = "fatal_mid_batch")
    )
  )
  S7::method(fit_model, FatalMidBatchFitter) <- function(
    fitter,
    data_bundle,
    fit_spec,
    seed,
    task_ctx
  ) {
    if (identical(task_ctx$rep_idx, 2L)) {
      stop(bayesim_config_error("fatal mid-batch fitter failure"))
    }
    fit_model(
      LinearRegressionFitter(n_draws = fitter@n_draws),
      data_bundle,
      fit_spec,
      seed,
      task_ctx
    )
  }

  path <- file.path(withr::local_tempdir(), "fatal-mid-batch")
  cfg <- simulation_config(
    data_grid = data.frame(n = 10L),
    fit_grid = data.frame(model = "lm"),
    data_generator = gen,
    fitter = FatalMidBatchFitter(n_draws = 10L),
    metrics = list(),
    n_replicates = 6L,
    seed = 21L,
    result_path = path,
    checkpoint_every = 3L
  )

  expect_error(
    run_simulation(cfg, resume = "never", progress = FALSE, verbose = FALSE),
    class = "bayesim_config_error"
  )

  # The successful siblings of the fatal task were persisted before the
  # re-raise; the fatal task itself stayed pending (it has no outcome).
  checkpoint <- read_checkpoint(path)
  expect_false(is.null(checkpoint))
  persisted <- vapply(
    checkpoint$task_outcomes,
    function(x) x$task_id,
    character(1)
  )
  expect_equal(
    sort(persisted),
    sort(c("d001_f001_r00001", "d001_f001_r00003"))
  )
  expect_true(all(vapply(
    checkpoint$task_outcomes,
    function(x) identical(x$status, "success"),
    logical(1)
  )))
  expect_false(
    "d001_f001_r00002" %in%
      checkpoint$task_grid$task_id[
        checkpoint$task_grid$status %in% c("success", "failed")
      ]
  )

  # The emergency checkpoint is a valid resume point.
  resume_state <- load_for_resume(path, cfg)
  expect_equal(nrow(resume_state$prior_results), 2L)
  expect_length(resume_state$prior_task_results, 2L)
  expect_equal(
    resume_state$task_grid$status[
      resume_state$task_grid$task_id == "d001_f001_r00002"
    ],
    "pending"
  )
})
