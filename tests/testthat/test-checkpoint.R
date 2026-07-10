# test-checkpoint.R
# Tests for checkpoint/resume functionality
# testthat 3e style with describe/it blocks

# Load the package
library(bayesim)

# =============================================================================
# Helper functions for test setup
# =============================================================================

#' Create a minimal valid checkpoint data structure for testing
#'
#' @param checkpoint_dir Path to checkpoint directory
#' @param checkpoint_id Checkpoint ID (integer, e.g., 1)
#' @param task_grid Data frame with task grid
#' @param results_df Data frame with task results
#' @param config_fingerprint Configuration fingerprint hash
#'
#' @return Invisibly returns the checkpoint directory path
#'
#' @keywords internal
create_test_checkpoint <- function(
  checkpoint_dir,
  checkpoint_id = 1L,
  task_grid = NULL,
  results_df = NULL,
  config_fingerprint = "test_fingerprint_abc123"
) {
  cp_name <- sprintf("cp_%06d", checkpoint_id)
  cp_path <- file.path(checkpoint_dir, "checkpoints", cp_name)
  dir.create(cp_path, recursive = TRUE, showWarnings = FALSE)

  # Create meta.json matching actual write_checkpoint format
  if (is.null(task_grid)) {
    task_grid <- data.frame(
      task_id = c("d001_f001_r00001", "d001_f001_r00002"),
      status = c("success", "pending"),
      stringsAsFactors = FALSE
    )
  }

  meta <- list(
    checkpoint_id = checkpoint_id,
    created = as.character(Sys.time()),
    config_fingerprint = config_fingerprint,
    n_tasks = nrow(task_grid),
    n_success = sum(task_grid$status == "success", na.rm = TRUE),
    n_failed = sum(task_grid$status == "failed", na.rm = TRUE),
    n_pending = sum(task_grid$status == "pending", na.rm = TRUE)
  )
  jsonlite::write_json(meta, file.path(cp_path, "meta.json"), auto_unbox = TRUE)

  # Create ledger.rds
  saveRDS(task_grid, file.path(cp_path, "ledger.rds"))

  # Create results.rds
  if (is.null(results_df)) {
    results_df <- data.frame(
      task_id = "d001_f001_r00001",
      status = "success",
      metric_rmse = 0.05,
      stringsAsFactors = FALSE
    )
  }
  saveRDS(results_df, file.path(cp_path, "results.rds"))

  # Create checksums.json
  checksums <- list(
    "meta.json" = digest::digest(file.path(cp_path, "meta.json"), file = TRUE),
    "ledger.rds" = digest::digest(
      file.path(cp_path, "ledger.rds"),
      file = TRUE
    ),
    "results.rds" = digest::digest(
      file.path(cp_path, "results.rds"),
      file = TRUE
    )
  )
  jsonlite::write_json(
    checksums,
    file.path(cp_path, "checksums.json"),
    auto_unbox = TRUE
  )

  invisible(checkpoint_dir)
}

#' Create a minimal run manifest for testing
#'
#' @param checkpoint_dir Path to checkpoint directory
#' @param config_fingerprint Configuration fingerprint hash
#'
#' @return Invisibly returns the checkpoint directory path
#'
#' @keywords internal
create_test_run_manifest <- function(
  checkpoint_dir,
  config_fingerprint = "test_fingerprint_abc123"
) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list(
    run_schema_version = 1L,
    result_schema_version = 1L,
    config_fingerprint = config_fingerprint,
    created = as.character(Sys.time())
  )
  jsonlite::write_json(
    manifest,
    file.path(checkpoint_dir, "run_manifest.json"),
    auto_unbox = TRUE
  )
  invisible(checkpoint_dir)
}

#' Create a latest.json pointer for testing
#'
#' @param checkpoint_dir Path to checkpoint directory
#' @param checkpoint_id Checkpoint ID to point to (integer or NULL)
#'
#' @return Invisibly returns the checkpoint directory path
#'
#' @keywords internal
create_test_latest_pointer <- function(checkpoint_dir, checkpoint_id = NULL) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  latest <- list(checkpoint_id = checkpoint_id)
  jsonlite::write_json(
    latest,
    file.path(checkpoint_dir, "latest.json"),
    auto_unbox = TRUE,
    null = "null"
  )
  invisible(checkpoint_dir)
}

#' Create a minimal task result object for testing
#'
#' @param task_id Task identifier
#' @param status Task status
#' @param metrics Named list of metrics
#'
#' @return A bayesim_task_result-like list
#'
#' @keywords internal
create_test_task_result <- function(
  task_id = "t1",
  status = "success",
  metrics = NULL
) {
  # Default metrics - use empty list for failed/skipped to avoid rbind column mismatch
  if (is.null(metrics)) {
    if (status == "success") {
      metrics <- list(rmse = 0.05)
    } else {
      metrics <- list()
    }
  }
  list(
    task_id = task_id,
    status = status,
    metrics = metrics,
    diagnostics = NULL,
    error = NULL,
    timing = list(total = 1.0)
  )
}

# =============================================================================
# 1. Checkpoint Directory Initialization
# =============================================================================

describe("Checkpoint Directory Initialization", {
  describe("init_checkpoint_dir()", {
    it("creates correct directory structure", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      expect_true(dir.exists(result_path))
      expect_true(dir.exists(file.path(result_path, "checkpoints")))
      expect_true(file.exists(file.path(result_path, "run_manifest.json")))
      expect_true(file.exists(file.path(result_path, "latest.json")))
    })

    it("run_manifest.json has correct schema version", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      manifest <- jsonlite::read_json(file.path(
        result_path,
        "run_manifest.json"
      ))
      expect_equal(manifest$run_schema_version, 1L)
      expect_equal(manifest$result_schema_version, 1L)
    })

    it("run_manifest.json includes config_fingerprint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(
        result_path,
        config_fingerprint = "my_custom_fingerprint"
      )

      manifest <- jsonlite::read_json(file.path(
        result_path,
        "run_manifest.json"
      ))
      expect_equal(manifest$config_fingerprint, "my_custom_fingerprint")
    })

    it("run_manifest.json includes created timestamp", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      manifest <- jsonlite::read_json(file.path(
        result_path,
        "run_manifest.json"
      ))
      expect_true(!is.null(manifest$created))
    })

    it("latest.json initialized with NULL checkpoint_id", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      latest <- jsonlite::read_json(file.path(result_path, "latest.json"))
      # jsonlite may return NULL or empty list for JSON null
      expect_true(
        is.null(latest$checkpoint_id) || length(latest$checkpoint_id) == 0
      )
    })

    it("creates checkpoints subdirectory", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      expect_true(dir.exists(file.path(result_path, "checkpoints")))
    })

    it("is idempotent - calling twice does not error", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      expect_silent(init_checkpoint_dir(
        result_path,
        config_fingerprint = "test_fingerprint"
      ))
      expect_silent(init_checkpoint_dir(
        result_path,
        config_fingerprint = "test_fingerprint"
      ))
    })

    it("returns NULL when result_path is NULL", {
      result <- init_checkpoint_dir(
        NULL,
        config_fingerprint = "test_fingerprint"
      )
      expect_null(result)
    })
  })
})

# =============================================================================
# 2. Checkpoint Writing
# =============================================================================

describe("Checkpoint Writing", {
  describe("get_next_checkpoint_id()", {
    it("returns 1 for empty checkpoints directory", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      dir.create(file.path(result_path, "checkpoints"), recursive = TRUE)

      result <- get_next_checkpoint_id(result_path)

      expect_equal(result, 1L)
    })

    it("returns correct sequence after one checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      dir.create(
        file.path(result_path, "checkpoints", "cp_000001"),
        recursive = TRUE
      )

      result <- get_next_checkpoint_id(result_path)

      expect_equal(result, 2L)
    })

    it("returns correct sequence after multiple checkpoints", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      dir.create(
        file.path(result_path, "checkpoints", "cp_000001"),
        recursive = TRUE
      )
      dir.create(
        file.path(result_path, "checkpoints", "cp_000002"),
        recursive = TRUE
      )
      dir.create(
        file.path(result_path, "checkpoints", "cp_000005"),
        recursive = TRUE
      )

      result <- get_next_checkpoint_id(result_path)

      # Should find highest existing and increment
      expect_equal(result, 6L)
    })

    it("ignores non-checkpoint directories", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      dir.create(
        file.path(result_path, "checkpoints", "cp_000001"),
        recursive = TRUE
      )
      dir.create(
        file.path(result_path, "checkpoints", "other_dir"),
        recursive = TRUE
      )
      dir.create(
        file.path(result_path, "checkpoints", "cp_000001.tmp"),
        recursive = TRUE
      )

      result <- get_next_checkpoint_id(result_path)

      expect_equal(result, 2L)
    })

    it("ignores .tmp directories when computing next ID", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      dir.create(
        file.path(result_path, "checkpoints", "cp_000001.tmp"),
        recursive = TRUE
      )

      result <- get_next_checkpoint_id(result_path)

      expect_equal(result, 1L)
    })
  })

  describe("write_checkpoint()", {
    it("prunes snapshots older than the configured retention count", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")
      task_grid <- data.frame(task_id = "t1", status = "success")
      task_results <- list(create_test_task_result())

      for (i in seq_len(4L)) {
        write_checkpoint(
          result_path,
          task_grid,
          task_results,
          config_fingerprint = "test_fingerprint",
          keep_checkpoints = 2L
        )
      }

      expect_equal(list_checkpoints(result_path), c(3L, 4L))
      expect_equal(
        jsonlite::read_json(file.path(
          result_path,
          "latest.json"
        ))$checkpoint_id,
        4L
      )
    })

    it("uses cached prior results without re-reading the previous checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")
      task_grid <- data.frame(task_id = "t1", status = "success")
      task_results <- list(create_test_task_result())
      local_mocked_bindings(
        read_checkpoint = function(...) stop("unexpected checkpoint read"),
        .package = "bayesim"
      )

      expect_no_error(write_checkpoint(
        result_path,
        task_grid,
        task_results,
        config_fingerprint = "test_fingerprint",
        prior_results_df = data.frame(
          task_id = character(),
          status = character()
        )
      ))
    })

    it("creates checkpoint directory with correct name", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "d001_f001_r00001",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(
        create_test_task_result(
          task_id = "d001_f001_r00001",
          status = "success"
        )
      )

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp"
      )

      expect_true(dir.exists(file.path(
        result_path,
        "checkpoints",
        "cp_000001"
      )))
    })

    it("creates meta.json with correct fields", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "pending"),
        stringsAsFactors = FALSE
      )
      task_results <- list(
        create_test_task_result(task_id = "t1", status = "success")
      )

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp_123"
      )

      meta <- jsonlite::read_json(file.path(
        result_path,
        "checkpoints",
        "cp_000001",
        "meta.json"
      ))
      expect_equal(meta$checkpoint_id, 1L)
      expect_equal(meta$config_fingerprint, "test_fp_123")
      expect_equal(meta$n_tasks, 2)
      expect_equal(meta$n_success, 1)
      expect_equal(meta$n_pending, 1)
    })

    it("creates ledger.rds with correct data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "pending"),
        stringsAsFactors = FALSE
      )
      task_results <- list(
        create_test_task_result(task_id = "t1", status = "success")
      )

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp"
      )

      saved_grid <- readRDS(file.path(
        result_path,
        "checkpoints",
        "cp_000001",
        "ledger.rds"
      ))
      expect_equal(saved_grid$task_id, c("t1", "t2"))
      expect_equal(saved_grid$status, c("success", "pending"))
    })

    it("creates results.rds with correct data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(
        create_test_task_result(
          task_id = "t1",
          status = "success",
          metrics = list(rmse = 0.05, bias = 0.01)
        )
      )

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp"
      )

      saved_results <- readRDS(file.path(
        result_path,
        "checkpoints",
        "cp_000001",
        "results.rds"
      ))
      expect_equal(saved_results$task_id, "t1")
      expect_equal(saved_results$rmse, 0.05)
      expect_equal(saved_results$bias, 0.01)
    })

    it("creates checksums.json with correct files", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(create_test_task_result())

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp"
      )

      checksums <- jsonlite::read_json(file.path(
        result_path,
        "checkpoints",
        "cp_000001",
        "checksums.json"
      ))
      expect_true("meta.json" %in% names(checksums))
      expect_true("ledger.rds" %in% names(checksums))
      expect_true("results.rds" %in% names(checksums))
    })

    it("updates latest.json to point to new checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(create_test_task_result())

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp"
      )

      latest <- jsonlite::read_json(file.path(result_path, "latest.json"))
      expect_equal(latest$checkpoint_id, 1L)
    })

    it("uses atomic write (tmp dir then rename)", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(create_test_task_result())

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp"
      )

      # After write, no .tmp directory should exist
      checkpoints_dir <- file.path(result_path, "checkpoints")
      tmp_dirs <- list.dirs(
        checkpoints_dir,
        full.names = FALSE,
        recursive = FALSE
      )
      tmp_dirs <- tmp_dirs[grepl("\\.tmp$", tmp_dirs)]
      expect_equal(length(tmp_dirs), 0)
    })

    it("returns NULL when result_path is NULL", {
      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(create_test_task_result())

      result <- write_checkpoint(
        result_path = NULL,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "test_fp"
      )

      expect_null(result)
    })
  })
})

# =============================================================================
# 3. Checkpoint Reading
# =============================================================================

describe("Checkpoint Reading", {
  describe("read_checkpoint()", {
    it("returns correct data structure for valid checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = c("d001_f001_r00001", "d001_f001_r00002"),
        status = c("success", "pending"),
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = "d001_f001_r00001",
        status = "success",
        metric_rmse = 0.05,
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df
      )
      create_test_latest_pointer(result_path, 1L)

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)

      expect_true(is.list(cp_data))
      expect_true(!is.null(cp_data$meta))
      expect_true(!is.null(cp_data$task_grid))
      expect_true(!is.null(cp_data$results_df))
    })

    it("returns NULL for non-existent checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      create_test_run_manifest(result_path)
      dir.create(file.path(result_path, "checkpoints"), showWarnings = FALSE)

      result <- read_checkpoint(result_path, checkpoint_id = 999L)

      expect_null(result)
    })

    it("returns NULL for checkpoint with missing files", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      create_test_run_manifest(result_path)

      # Create incomplete checkpoint (missing ledger.rds)
      cp_path <- file.path(result_path, "checkpoints", "cp_000001")
      dir.create(cp_path, recursive = TRUE)
      jsonlite::write_json(
        list(checkpoint_id = 1L),
        file.path(cp_path, "meta.json"),
        auto_unbox = TRUE
      )

      result <- NULL
      expect_warning(
        result <- read_checkpoint(result_path, checkpoint_id = 1L),
        "invalid checksums"
      )

      expect_null(result)
    })

    it("validates checksums and succeeds on valid checksums", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df
      )

      # Should not throw or warn
      expect_silent(read_checkpoint(result_path, checkpoint_id = 1L))
    })

    it("warns on invalid checksums", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df
      )

      # Corrupt a checksum
      cp_path <- file.path(result_path, "checkpoints", "cp_000001")
      checksums <- jsonlite::read_json(file.path(cp_path, "checksums.json"))
      checksums["meta.json"] <- "invalid_hash"
      jsonlite::write_json(
        checksums,
        file.path(cp_path, "checksums.json"),
        auto_unbox = TRUE
      )

      expect_warning(
        read_checkpoint(result_path, checkpoint_id = 1L),
        "checksum"
      )
    })

    it("returns correct meta information", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df,
        config_fingerprint = "my_special_fingerprint"
      )

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)

      expect_equal(cp_data$meta$checkpoint_id, 1L)
      expect_equal(cp_data$meta$config_fingerprint, "my_special_fingerprint")
    })

    it("returns correct task_grid data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = c("t1", "t2", "t3"),
        status = c("success", "failed", "pending"),
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df
      )

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)

      expect_equal(nrow(cp_data$task_grid), 3)
      expect_equal(cp_data$task_grid$status, c("success", "failed", "pending"))
    })

    it("returns correct results_df data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "pending"),
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = c("t1", "t2"),
        metric_rmse = c(0.05, 0.10),
        metric_bias = c(0.01, 0.02),
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df
      )

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)

      expect_equal(nrow(cp_data$results_df), 2)
      expect_equal(cp_data$results_df$metric_rmse, c(0.05, 0.10))
    })

    it("reads latest checkpoint when checkpoint_id is NULL", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df
      )
      create_test_latest_pointer(result_path, 1L)

      cp_data <- read_checkpoint(result_path)

      expect_equal(cp_data$checkpoint_id, 1L)
    })

    it("returns NULL when result_path is NULL", {
      result <- read_checkpoint(NULL)
      expect_null(result)
    })
  })
})

# =============================================================================
# 4. Checkpoint Validation
# =============================================================================

describe("Checkpoint Validation", {
  describe("validate_checkpoint_fingerprint()", {
    it("returns TRUE for matching fingerprints", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(
        result_path,
        config_fingerprint = "test_fp_match"
      )
      create_test_checkpoint(result_path, config_fingerprint = "test_fp_match")

      checkpoint <- read_checkpoint(result_path, checkpoint_id = 1L)
      result <- validate_checkpoint_fingerprint(checkpoint, "test_fp_match")

      expect_true(result)
    })

    it("returns FALSE for mismatched fingerprints", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(
        result_path,
        config_fingerprint = "test_fp_original"
      )
      create_test_checkpoint(
        result_path,
        config_fingerprint = "test_fp_original"
      )

      checkpoint <- read_checkpoint(result_path, checkpoint_id = 1L)
      result <- validate_checkpoint_fingerprint(
        checkpoint,
        "different_fingerprint"
      )

      expect_false(result)
    })

    it("returns FALSE for NULL checkpoint", {
      result <- validate_checkpoint_fingerprint(NULL, "any_fingerprint")
      expect_false(result)
    })

    it("returns FALSE for checkpoint missing meta", {
      checkpoint <- list(checkpoint_id = 1L) # no meta
      result <- validate_checkpoint_fingerprint(checkpoint, "any_fingerprint")
      expect_false(result)
    })
  })
})

# =============================================================================
# 5. Resume Logic
# =============================================================================

describe("Resume Logic", {
  describe("can_resume()", {
    it("returns TRUE for valid run with checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      results_df <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        task_grid = task_grid,
        results_df = results_df
      )
      create_test_latest_pointer(result_path, 1L)

      expect_true(can_resume(result_path))
    })

    it("returns FALSE for missing manifest", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      dir.create(result_path, showWarnings = FALSE)

      expect_false(can_resume(result_path))
    })

    it("returns FALSE for empty results directory", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      expect_false(can_resume(result_path))
    })

    it("returns FALSE for run with no checkpoints", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      dir.create(file.path(result_path, "checkpoints"), showWarnings = FALSE)
      create_test_latest_pointer(result_path, NULL)

      expect_false(can_resume(result_path))
    })

    it("returns FALSE for corrupted latest.json", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      dir.create(file.path(result_path, "checkpoints"), showWarnings = FALSE)

      # Write invalid JSON
      writeLines("not valid json", file.path(result_path, "latest.json"))

      expect_false(can_resume(result_path))
    })

    it("returns FALSE when result_path is NULL", {
      expect_false(can_resume(NULL))
    })
  })

  describe("get_latest_valid_checkpoint()", {
    it("finds newest valid checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)

      # Create multiple checkpoints
      for (i in 1:3) {
        create_test_checkpoint(result_path, checkpoint_id = i)
      }
      create_test_latest_pointer(result_path, 3L)

      result <- get_latest_valid_checkpoint(result_path)

      expect_true(is.list(result))
      expect_equal(result$checkpoint_id, 3L)
    })

    it("falls back to previous checkpoint if latest is invalid", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)

      # Create two valid checkpoints
      create_test_checkpoint(result_path, checkpoint_id = 1L)
      create_test_checkpoint(result_path, checkpoint_id = 2L)

      # Create invalid checkpoint (missing required files)
      cp_path <- file.path(result_path, "checkpoints", "cp_000003")
      dir.create(cp_path, recursive = TRUE)
      # No meta.json, ledger.rds, or results.rds

      create_test_latest_pointer(result_path, 3L)

      result <- NULL
      expect_warning(
        result <- get_latest_valid_checkpoint(result_path),
        "invalid checksums"
      )

      expect_equal(result$checkpoint_id, 2L)
    })

    it("returns NULL when no valid checkpoints exist", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      dir.create(file.path(result_path, "checkpoints"), showWarnings = FALSE)
      create_test_latest_pointer(result_path, 1L)

      # No actual checkpoint files

      result <- get_latest_valid_checkpoint(result_path)

      expect_null(result)
    })

    it("ignores .tmp directories when finding valid checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, checkpoint_id = 1L)

      # Create .tmp directory (incomplete checkpoint)
      tmp_path <- file.path(result_path, "checkpoints", "cp_000002.tmp")
      dir.create(tmp_path, recursive = TRUE)

      create_test_latest_pointer(result_path, 1L)

      result <- get_latest_valid_checkpoint(result_path)

      expect_equal(result$checkpoint_id, 1L)
    })

    it("validates fingerprint when provided", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        checkpoint_id = 1L,
        config_fingerprint = "matching_fp"
      )
      create_test_latest_pointer(result_path, 1L)

      result <- get_latest_valid_checkpoint(
        result_path,
        config_fingerprint = "matching_fp"
      )

      expect_equal(result$checkpoint_id, 1L)
    })

    it("returns NULL when fingerprint does not match", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        checkpoint_id = 1L,
        config_fingerprint = "original_fp"
      )
      create_test_latest_pointer(result_path, 1L)

      result <- get_latest_valid_checkpoint(
        result_path,
        config_fingerprint = "different_fp"
      )

      expect_null(result)
    })

    it("returns NULL when result_path is NULL", {
      result <- get_latest_valid_checkpoint(NULL)
      expect_null(result)
    })
  })

  describe("merge_task_grid_status()", {
    it("merges status from checkpoint into task grid", {
      # Original task grid
      task_grid <- data.frame(
        task_id = c("t1", "t2", "t3", "t4"),
        status = c("pending", "pending", "pending", "pending"),
        stringsAsFactors = FALSE
      )

      # Checkpoint ledger with some completed tasks
      checkpoint_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "failed"),
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, checkpoint_grid)

      expect_equal(result$status, c("success", "failed", "pending", "pending"))
    })

    it("preserves task grid columns not in checkpoint grid", {
      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("pending", "pending"),
        data_idx = c(1L, 1L),
        fit_idx = c(1L, 2L),
        stringsAsFactors = FALSE
      )

      checkpoint_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, checkpoint_grid)

      expect_equal(result$data_idx, c(1L, 1L))
      expect_equal(result$fit_idx, c(1L, 2L))
    })

    it("handles empty checkpoint grid (all pending)", {
      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("pending", "pending"),
        stringsAsFactors = FALSE
      )

      checkpoint_grid <- data.frame(
        task_id = character(),
        status = character(),
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, checkpoint_grid)

      expect_equal(result$status, c("pending", "pending"))
    })

    it("handles task not in grid (skip silently)", {
      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("pending", "pending"),
        stringsAsFactors = FALSE
      )

      # Checkpoint grid has task that doesn't exist in fresh grid
      checkpoint_grid <- data.frame(
        task_id = c("t1", "t999"),
        status = c("success", "success"),
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, checkpoint_grid)

      expect_equal(result$status, c("success", "pending"))
    })
  })

  describe("merge_results()", {
    it("accepts byte-identical duplicate rows", {
      old_results <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "success"),
        metric_rmse = c(0.05, 0.10),
        stringsAsFactors = FALSE
      )

      new_results <- data.frame(
        task_id = c("t2", "t3"),
        status = c("success", "success"),
        metric_rmse = c(0.10, 0.15),
        stringsAsFactors = FALSE
      )

      result <- merge_results(old_results, new_results)

      expect_equal(nrow(result), 3)
      expect_equal(result$metric_rmse[result$task_id == "t2"], 0.10)
    })

    it("errors on conflicting duplicate terminal rows", {
      old_results <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "success"),
        metric_rmse = c(0.05, 0.10),
        stringsAsFactors = FALSE
      )

      new_results <- data.frame(
        task_id = c("t2", "t3"),
        status = c("success", "success"),
        metric_rmse = c(0.08, 0.15),
        stringsAsFactors = FALSE
      )

      expect_error(
        merge_results(old_results, new_results),
        "Conflicting duplicate terminal rows"
      )
    })

    it("preserves unique tasks from both sources", {
      old_results <- data.frame(
        task_id = c("t1", "t2"),
        metric_rmse = c(0.05, 0.10),
        stringsAsFactors = FALSE
      )

      new_results <- data.frame(
        task_id = c("t3", "t4"),
        metric_rmse = c(0.15, 0.20),
        stringsAsFactors = FALSE
      )

      result <- merge_results(old_results, new_results)

      expect_equal(nrow(result), 4)
      expect_equal(sort(result$task_id), c("t1", "t2", "t3", "t4"))
    })

    it("handles empty old results", {
      old_results <- data.frame(
        task_id = character(),
        metric_rmse = numeric(),
        stringsAsFactors = FALSE
      )

      new_results <- data.frame(
        task_id = "t1",
        metric_rmse = 0.05,
        stringsAsFactors = FALSE
      )

      result <- merge_results(old_results, new_results)

      expect_equal(nrow(result), 1)
      expect_equal(result$task_id, "t1")
    })

    it("handles NULL old results", {
      new_results <- data.frame(
        task_id = "t1",
        metric_rmse = 0.05,
        stringsAsFactors = FALSE
      )

      result <- merge_results(NULL, new_results)

      expect_equal(nrow(result), 1)
      expect_equal(result$task_id, "t1")
    })

    it("handles NULL new results", {
      old_results <- data.frame(
        task_id = "t1",
        metric_rmse = 0.05,
        stringsAsFactors = FALSE
      )

      result <- merge_results(old_results, NULL)

      expect_equal(nrow(result), 1)
      expect_equal(result$task_id, "t1")
    })

    it("handles both NULL results", {
      result <- merge_results(NULL, NULL)
      expect_null(result)
    })
  })
})

# =============================================================================
# 6. Integration Tests
# =============================================================================

describe("Checkpoint Integration", {
  describe("Full write -> read cycle", {
    it("preserves data through write and read", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      # Original data
      original_task_grid <- data.frame(
        task_id = c("d001_f001_r00001", "d001_f001_r00002", "d001_f001_r00003"),
        status = c("success", "failed", "pending"),
        stringsAsFactors = FALSE
      )
      # Use consistent metrics structure for all results to avoid rbind column mismatch
      original_task_results <- list(
        create_test_task_result(
          task_id = "d001_f001_r00001",
          status = "success",
          metrics = list(rmse = 0.05, bias = 0.01)
        ),
        create_test_task_result(
          task_id = "d001_f001_r00002",
          status = "failed",
          metrics = list(rmse = NA_real_, bias = NA_real_)
        )
      )

      # Write checkpoint
      write_checkpoint(
        result_path = result_path,
        task_grid = original_task_grid,
        task_results = original_task_results,
        config_fingerprint = "integration_test_fp"
      )

      # Read checkpoint
      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)

      # Verify data preserved
      expect_equal(cp_data$task_grid$task_id, original_task_grid$task_id)
      expect_equal(cp_data$task_grid$status, original_task_grid$status)
      expect_true("d001_f001_r00001" %in% cp_data$results_df$task_id)
    })

    it("preserves numeric precision", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      # Data with specific numeric values
      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(
        create_test_task_result(
          task_id = "t1",
          status = "success",
          metrics = list(
            rmse = 0.123456789,
            bias = -0.987654321,
            coverage = 0.95
          )
        )
      )

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "precision_test"
      )

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)

      expect_equal(cp_data$results_df$rmse, 0.123456789, tolerance = 1e-10)
      expect_equal(cp_data$results_df$bias, -0.987654321, tolerance = 1e-10)
    })
  })

  describe("Multiple checkpoints in sequence", {
    it("creates sequential checkpoint IDs", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(create_test_task_result())

      # Write three checkpoints
      for (i in 1:3) {
        write_checkpoint(
          result_path = result_path,
          task_grid = task_grid,
          task_results = task_results,
          config_fingerprint = "seq_test"
        )
      }

      # Verify all three exist
      expect_true(dir.exists(file.path(
        result_path,
        "checkpoints",
        "cp_000001"
      )))
      expect_true(dir.exists(file.path(
        result_path,
        "checkpoints",
        "cp_000002"
      )))
      expect_true(dir.exists(file.path(
        result_path,
        "checkpoints",
        "cp_000003"
      )))

      # Verify latest.json points to newest
      latest <- jsonlite::read_json(file.path(result_path, "latest.json"))
      expect_equal(latest$checkpoint_id, 3L)
    })

    it("can read any checkpoint in sequence", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      # Create checkpoints with different data
      for (i in 1:3) {
        task_grid <- data.frame(
          task_id = "t1",
          status = "success",
          stringsAsFactors = FALSE
        )
        task_results <- list(
          create_test_task_result(
            task_id = "t1",
            status = "success",
            metrics = list(value = i)
          )
        )

        write_checkpoint(
          result_path = result_path,
          task_grid = task_grid,
          task_results = task_results,
          config_fingerprint = "seq_read_test"
        )
      }

      # Read each and verify correct data
      cp1 <- read_checkpoint(result_path, checkpoint_id = 1L)
      cp2 <- read_checkpoint(result_path, checkpoint_id = 2L)
      cp3 <- read_checkpoint(result_path, checkpoint_id = 3L)

      expect_equal(cp1$results_df$value, 1)
      expect_equal(cp2$results_df$value, 2)
      expect_equal(cp3$results_df$value, 3)
    })

    it("get_latest_valid_checkpoint finds newest after multiple writes", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(create_test_task_result())

      for (i in 1:5) {
        write_checkpoint(
          result_path = result_path,
          task_grid = task_grid,
          task_results = task_results,
          config_fingerprint = "find_test"
        )
      }

      result <- get_latest_valid_checkpoint(result_path)

      expect_equal(result$checkpoint_id, 5L)
    })

    it("can resume from intermediate checkpoint if latest is corrupted", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      # Create two valid checkpoints
      for (i in 1:2) {
        task_grid <- data.frame(
          task_id = "t1",
          status = "success",
          stringsAsFactors = FALSE
        )
        task_results <- list(
          create_test_task_result(
            task_id = "t1",
            status = "success",
            metrics = list(value = i)
          )
        )

        write_checkpoint(
          result_path = result_path,
          task_grid = task_grid,
          task_results = task_results,
          config_fingerprint = "corruption_test"
        )
      }

      # Corrupt the latest checkpoint
      meta_path <- file.path(
        result_path,
        "checkpoints",
        "cp_000002",
        "meta.json"
      )
      file.remove(meta_path)

      # get_latest_valid_checkpoint should return the previous valid one
      valid_cp <- NULL
      expect_warning(
        valid_cp <- get_latest_valid_checkpoint(result_path),
        "invalid checksums"
      )
      expect_equal(valid_cp$checkpoint_id, 1L)

      # And we should be able to read it
      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)
      expect_equal(cp_data$results_df$value, 1)
    })
  })

  describe("list_checkpoints()", {
    it("returns empty integer vector for no checkpoints", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      result <- list_checkpoints(result_path)

      expect_equal(result, integer(0))
    })

    it("returns sorted checkpoint IDs", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(create_test_task_result())

      for (i in c(1, 3, 5)) {
        # Create checkpoint manually to test non-sequential IDs
        create_test_checkpoint(result_path, checkpoint_id = i)
      }

      result <- list_checkpoints(result_path)

      expect_equal(result, c(1L, 3L, 5L))
    })
  })

  describe("Edge cases", {
    it("handles empty task_grid", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = character(),
        status = character(),
        stringsAsFactors = FALSE
      )
      task_results <- list()

      expect_silent(
        write_checkpoint(
          result_path = result_path,
          task_grid = task_grid,
          task_results = task_results,
          config_fingerprint = "empty_test"
        )
      )

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)
      expect_equal(nrow(cp_data$task_grid), 0)
      expect_equal(nrow(cp_data$results_df), 0)
    })

    it("handles large task IDs", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      # Very large task ID
      task_grid <- data.frame(
        task_id = "d999_f999_r99999",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(
        create_test_task_result(
          task_id = "d999_f999_r99999",
          status = "success"
        )
      )

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "large_id_test"
      )

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)
      expect_equal(cp_data$task_grid$task_id, "d999_f999_r99999")
    })

    it("handles special characters in metric names", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path, config_fingerprint = "test_fingerprint")

      task_grid <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )
      task_results <- list(
        create_test_task_result(
          task_id = "t1",
          status = "success",
          metrics = list("metric_with_special" = 1.0)
        )
      )

      write_checkpoint(
        result_path = result_path,
        task_grid = task_grid,
        task_results = task_results,
        config_fingerprint = "special_chars_test"
      )

      cp_data <- read_checkpoint(result_path, checkpoint_id = 1L)
      expect_true("metric_with_special" %in% names(cp_data$results_df))
    })
  })
})
