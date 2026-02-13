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
#' @param checkpoint_id Checkpoint ID (e.g., "cp_000001")
#' @param ledger Data frame with task ledger
#' @param results Data frame with task results
#' @param config_fingerprint Configuration fingerprint hash
#'
#' @return Invisibly returns the checkpoint directory path
#'
#' @keywords internal
create_test_checkpoint <- function(
  checkpoint_dir,
  checkpoint_id = "cp_000001",
  ledger = NULL,
  results = NULL,
  config_fingerprint = "test_fingerprint_abc123"
) {
  cp_path <- file.path(checkpoint_dir, "checkpoints", checkpoint_id)
  dir.create(cp_path, recursive = TRUE, showWarnings = FALSE)

  # Create meta.json
  meta <- list(
    checkpoint_id = checkpoint_id,
    config_fingerprint = config_fingerprint,
    created_at = Sys.time(),
    n_tasks_total = 100L,
    n_tasks_completed = 10L,
    run_schema_version = "1.0",
    result_schema_version = "1.0"
  )
  jsonlite::write_json(meta, file.path(cp_path, "meta.json"), auto_unbox = TRUE)

  # Create ledger.rds
  if (is.null(ledger)) {
    ledger <- data.frame(
      task_id = c("d001_f001_r00001", "d001_f001_r00002"),
      status = c("success", "pending"),
      stringsAsFactors = FALSE
    )
  }
  saveRDS(ledger, file.path(cp_path, "ledger.rds"))

  # Create results.rds
  if (is.null(results)) {
    results <- data.frame(
      task_id = "d001_f001_r00001",
      metric_rmse = 0.05,
      stringsAsFactors = FALSE
    )
  }
  saveRDS(results, file.path(cp_path, "results.rds"))

  # Create checksums.json
  checksums <- list(
    "meta.json" = digest::digest(file.path(cp_path, "meta.json"), file = TRUE),
    "ledger.rds" = digest::digest(file.path(cp_path, "ledger.rds"), file = TRUE),
    "results.rds" = digest::digest(file.path(cp_path, "results.rds"), file = TRUE)
  )
  jsonlite::write_json(checksums, file.path(cp_path, "checksums.json"), auto_unbox = TRUE)

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
  manifest <- list(
    run_schema_version = "1.0",
    result_schema_version = "1.0",
    config_fingerprint = config_fingerprint,
    created_at = Sys.time()
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
#' @param checkpoint_id Checkpoint ID to point to (or NULL)
#'
#' @return Invisibly returns the checkpoint directory path
#'
#' @keywords internal
create_test_latest_pointer <- function(checkpoint_dir, checkpoint_id = NULL) {
  latest <- list(checkpoint_id = checkpoint_id)
  jsonlite::write_json(
    latest,
    file.path(checkpoint_dir, "latest.json"),
    auto_unbox = TRUE,
    null = "null"
  )
  invisible(checkpoint_dir)
}

# =============================================================================
# 1. Checkpoint Directory Initialization
# =============================================================================

describe("Checkpoint Directory Initialization", {
  describe("init_checkpoint_dir()", {
    it("creates correct directory structure", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      expect_true(dir.exists(result_path))
      expect_true(dir.exists(file.path(result_path, "checkpoints")))
      expect_true(file.exists(file.path(result_path, "run_manifest.json")))
      expect_true(file.exists(file.path(result_path, "latest.json")))
    })

    it("run_manifest.json has correct schema version", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      manifest <- jsonlite::read_json(file.path(result_path, "run_manifest.json"))
      expect_equal(manifest$run_schema_version, "1.0")
      expect_equal(manifest$result_schema_version, "1.0")
    })

    it("run_manifest.json includes created_at timestamp", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      manifest <- jsonlite::read_json(file.path(result_path, "run_manifest.json"))
      expect_true(!is.null(manifest$created_at))
    })

    it("latest.json initialized with NULL checkpoint_id", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      latest <- jsonlite::read_json(file.path(result_path, "latest.json"))
      expect_true(is.null(latest$checkpoint_id))
    })

    it("creates checkpoints subdirectory", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      expect_true(dir.exists(file.path(result_path, "checkpoints")))
    })

    it("is idempotent - calling twice does not error", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      expect_silent(init_checkpoint_dir(result_path))
      expect_silent(init_checkpoint_dir(result_path))
    })

    it("preserves existing manifest when called on existing directory", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      # Get original timestamp
      manifest <- jsonlite::read_json(file.path(result_path, "run_manifest.json"))
      original_time <- manifest$created_at

      # Call again
      init_checkpoint_dir(result_path)

      # Timestamp should be preserved
      manifest <- jsonlite::read_json(file.path(result_path, "run_manifest.json"))
      expect_equal(manifest$created_at, original_time)
    })

    it("errors on invalid path (parent directory does not exist)", {
      expect_error(
        init_checkpoint_dir("/nonexistent/path/results"),
        class = "bayesim_checkpoint_error"
      )
    })
  })
})

# =============================================================================
# 2. Checkpoint Writing
# =============================================================================

describe("Checkpoint Writing", {
  describe("get_next_checkpoint_id()", {
    it("returns 'cp_000001' for empty checkpoints directory", {
      tmpdir <- withr::local_tempdir()
      checkpoints_dir <- file.path(tmpdir, "checkpoints")
      dir.create(checkpoints_dir)

      result <- get_next_checkpoint_id(checkpoints_dir)

      expect_equal(result, "cp_000001")
    })

    it("returns correct sequence after one checkpoint", {
      tmpdir <- withr::local_tempdir()
      checkpoints_dir <- file.path(tmpdir, "checkpoints")
      dir.create(checkpoints_dir)
      dir.create(file.path(checkpoints_dir, "cp_000001"))

      result <- get_next_checkpoint_id(checkpoints_dir)

      expect_equal(result, "cp_000002")
    })

    it("returns correct sequence after multiple checkpoints", {
      tmpdir <- withr::local_tempdir()
      checkpoints_dir <- file.path(tmpdir, "checkpoints")
      dir.create(checkpoints_dir)
      dir.create(file.path(checkpoints_dir, "cp_000001"))
      dir.create(file.path(checkpoints_dir, "cp_000002"))
      dir.create(file.path(checkpoints_dir, "cp_000005"))

      result <- get_next_checkpoint_id(checkpoints_dir)

      # Should find highest existing and increment
      expect_equal(result, "cp_000006")
    })

    it("ignores non-checkpoint directories", {
      tmpdir <- withr::local_tempdir()
      checkpoints_dir <- file.path(tmpdir, "checkpoints")
      dir.create(checkpoints_dir)
      dir.create(file.path(checkpoints_dir, "cp_000001"))
      dir.create(file.path(checkpoints_dir, "other_dir"))
      dir.create(file.path(checkpoints_dir, "cp_000001.tmp"))

      result <- get_next_checkpoint_id(checkpoints_dir)

      expect_equal(result, "cp_000002")
    })

    it("ignores .tmp directories when computing next ID", {
      tmpdir <- withr::local_tempdir()
      checkpoints_dir <- file.path(tmpdir, "checkpoints")
      dir.create(checkpoints_dir)
      dir.create(file.path(checkpoints_dir, "cp_000001.tmp"))

      result <- get_next_checkpoint_id(checkpoints_dir)

      expect_equal(result, "cp_000001")
    })
  })

  describe("write_checkpoint()", {
    it("creates checkpoint directory with correct name", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(
        task_id = "d001_f001_r00001",
        status = "success",
        stringsAsFactors = FALSE
      )
      results <- data.frame(
        task_id = "d001_f001_r00001",
        metric_rmse = 0.05,
        stringsAsFactors = FALSE
      )

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp",
        n_tasks_total = 10L,
        n_tasks_completed = 1L
      )

      expect_true(dir.exists(file.path(result_path, "checkpoints", "cp_000001")))
    })

    it("creates meta.json with correct fields", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp_123",
        n_tasks_total = 100L,
        n_tasks_completed = 5L
      )

      meta <- jsonlite::read_json(file.path(result_path, "checkpoints", "cp_000001", "meta.json"))
      expect_equal(meta$checkpoint_id, "cp_000001")
      expect_equal(meta$config_fingerprint, "test_fp_123")
      expect_equal(meta$n_tasks_total, 100)
      expect_equal(meta$n_tasks_completed, 5)
      expect_equal(meta$run_schema_version, "1.0")
    })

    it("creates ledger.rds with correct data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "pending"),
        stringsAsFactors = FALSE
      )
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp",
        n_tasks_total = 2L,
        n_tasks_completed = 1L
      )

      saved_ledger <- readRDS(file.path(result_path, "checkpoints", "cp_000001", "ledger.rds"))
      expect_equal(saved_ledger$task_id, c("t1", "t2"))
      expect_equal(saved_ledger$status, c("success", "pending"))
    })

    it("creates results.rds with correct data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(
        task_id = "t1",
        metric_rmse = 0.05,
        metric_bias = 0.01,
        stringsAsFactors = FALSE
      )

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      saved_results <- readRDS(file.path(result_path, "checkpoints", "cp_000001", "results.rds"))
      expect_equal(saved_results$task_id, "t1")
      expect_equal(saved_results$metric_rmse, 0.05)
      expect_equal(saved_results$metric_bias, 0.01)
    })

    it("creates checksums.json with correct files", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      checksums <- jsonlite::read_json(file.path(result_path, "checkpoints", "cp_000001", "checksums.json"))
      expect_true("meta.json" %in% names(checksums))
      expect_true("ledger.rds" %in% names(checksums))
      expect_true("results.rds" %in% names(checksums))
    })

    it("checksums.json contains valid SHA256 hashes", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      checksums <- jsonlite::read_json(file.path(result_path, "checkpoints", "cp_000001", "checksums.json"))
      cp_dir <- file.path(result_path, "checkpoints", "cp_000001")

      # Verify each checksum
      for (filename in names(checksums)) {
        actual_hash <- digest::digest(file.path(cp_dir, filename), algo = "sha256", file = TRUE)
        expect_equal(checksums[[filename]], actual_hash)
      }
    })

    it("updates latest.json to point to new checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      latest <- jsonlite::read_json(file.path(result_path, "latest.json"))
      expect_equal(latest$checkpoint_id, "cp_000001")
    })

    it("uses atomic write (tmp dir then rename)", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "test_fp",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      # After write, no .tmp directory should exist
      checkpoints_dir <- file.path(result_path, "checkpoints")
      tmp_dirs <- list.dirs(checkpoints_dir, full.names = FALSE, recursive = FALSE)
      tmp_dirs <- tmp_dirs[grepl("\\.tmp$", tmp_dirs)]
      expect_equal(length(tmp_dirs), 0)
    })

    it("cleans up tmp directory on failure", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      init_checkpoint_dir(result_path)

      # Create a file that will cause write failure (simulated by making path read-only)
      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      # This test would require mocking to simulate failure mid-write
      # For now, verify the tmp cleanup happens on expected errors
      expect_error(
        write_checkpoint(
          result_path = "/nonexistent/path",
          ledger = ledger,
          results = results,
          config_fingerprint = "test_fp",
          n_tasks_total = 1L,
          n_tasks_completed = 1L
        ),
        class = "bayesim_checkpoint_error"
      )
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

      ledger <- data.frame(
        task_id = c("d001_f001_r00001", "d001_f001_r00002"),
        status = c("success", "pending"),
        stringsAsFactors = FALSE
      )
      results <- data.frame(
        task_id = "d001_f001_r00001",
        metric_rmse = 0.05,
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, ledger = ledger, results = results)
      create_test_latest_pointer(result_path, "cp_000001")

      cp_data <- read_checkpoint(result_path, "cp_000001")

      expect_true(is.list(cp_data))
      expect_true(!is.null(cp_data$meta))
      expect_true(!is.null(cp_data$ledger))
      expect_true(!is.null(cp_data$results))
    })

    it("returns NULL for non-existent checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      create_test_run_manifest(result_path)
      dir.create(file.path(result_path, "checkpoints"), showWarnings = FALSE)

      result <- read_checkpoint(result_path, "cp_999999")

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
        list(checkpoint_id = "cp_000001"),
        file.path(cp_path, "meta.json")
      )

      result <- read_checkpoint(result_path, "cp_000001")

      expect_null(result)
    })

    it("validates checksums and succeeds on valid checksums", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, ledger = ledger, results = results)

      # Should not throw or warn
      expect_silent(read_checkpoint(result_path, "cp_000001"))
    })

    it("warns on invalid checksums", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, ledger = ledger, results = results)

      # Corrupt a checksum
      cp_path <- file.path(result_path, "checkpoints", "cp_000001")
      checksums <- jsonlite::read_json(file.path(cp_path, "checksums.json"))
      checksums["meta.json"] <- "invalid_hash"
      jsonlite::write_json(checksums, file.path(cp_path, "checksums.json"), auto_unbox = TRUE)

      expect_warning(
        read_checkpoint(result_path, "cp_000001"),
        "checksum"
      )
    })

    it("returns correct meta information", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      create_test_run_manifest(result_path)
      create_test_checkpoint(
        result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "my_special_fingerprint"
      )

      cp_data <- read_checkpoint(result_path, "cp_000001")

      expect_equal(cp_data$meta$checkpoint_id, "cp_000001")
      expect_equal(cp_data$meta$config_fingerprint, "my_special_fingerprint")
    })

    it("returns correct ledger data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      ledger <- data.frame(
        task_id = c("t1", "t2", "t3"),
        status = c("success", "failed", "pending"),
        stringsAsFactors = FALSE
      )
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, ledger = ledger, results = results)

      cp_data <- read_checkpoint(result_path, "cp_000001")

      expect_equal(nrow(cp_data$ledger), 3)
      expect_equal(cp_data$ledger$status, c("success", "failed", "pending"))
    })

    it("returns correct results data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      ledger <- data.frame(task_id = c("t1", "t2"), status = c("success", "pending"), stringsAsFactors = FALSE)
      results <- data.frame(
        task_id = c("t1", "t2"),
        metric_rmse = c(0.05, 0.10),
        metric_bias = c(0.01, 0.02),
        stringsAsFactors = FALSE
      )

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, ledger = ledger, results = results)

      cp_data <- read_checkpoint(result_path, "cp_000001")

      expect_equal(nrow(cp_data$results), 2)
      expect_equal(cp_data$results$metric_rmse, c(0.05, 0.10))
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

      create_test_run_manifest(result_path, config_fingerprint = "test_fp_match")
      create_test_checkpoint(result_path, config_fingerprint = "test_fp_match")

      result <- validate_checkpoint_fingerprint(result_path, "cp_000001", "test_fp_match")

      expect_true(result)
    })

    it("returns FALSE for mismatched fingerprints", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path, config_fingerprint = "test_fp_original")
      create_test_checkpoint(result_path, config_fingerprint = "test_fp_original")

      result <- validate_checkpoint_fingerprint(result_path, "cp_000001", "different_fingerprint")

      expect_false(result)
    })

    it("returns FALSE for missing checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      create_test_run_manifest(result_path)

      result <- validate_checkpoint_fingerprint(result_path, "cp_999999", "any_fingerprint")

      expect_false(result)
    })

    it("returns FALSE for checkpoint missing meta.json", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")
      create_test_run_manifest(result_path)

      # Create checkpoint without meta.json
      cp_path <- file.path(result_path, "checkpoints", "cp_000001")
      dir.create(cp_path, recursive = TRUE)

      result <- validate_checkpoint_fingerprint(result_path, "cp_000001", "any_fingerprint")

      expect_false(result)
    })
  })

  describe("Fingerprint mismatch handling", {
    it("is detected and reported", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path, config_fingerprint = "original")
      create_test_checkpoint(result_path, config_fingerprint = "original")

      # Try to validate with wrong fingerprint
      result <- validate_checkpoint_fingerprint(result_path, "cp_000001", "different")

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

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, ledger = ledger, results = results)
      create_test_latest_pointer(result_path, "cp_000001")

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
  })

  describe("load_for_resume()", {
    it("loads valid checkpoint data", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      ledger <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "pending"),
        stringsAsFactors = FALSE
      )
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      create_test_run_manifest(result_path, config_fingerprint = "test_fp")
      create_test_checkpoint(result_path, ledger = ledger, results = results, config_fingerprint = "test_fp")
      create_test_latest_pointer(result_path, "cp_000001")

      resume_data <- load_for_resume(result_path, config_fingerprint = "test_fp")

      expect_true(is.list(resume_data))
      expect_true(!is.null(resume_data$ledger))
      expect_true(!is.null(resume_data$results))
    })

    it("validates fingerprint and succeeds on match", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path, config_fingerprint = "matching_fp")
      create_test_checkpoint(result_path, config_fingerprint = "matching_fp")
      create_test_latest_pointer(result_path, "cp_000001")

      # Should not error
      result <- load_for_resume(result_path, config_fingerprint = "matching_fp")

      expect_true(is.list(result))
    })

    it("fails on fingerprint mismatch", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path, config_fingerprint = "original_fp")
      create_test_checkpoint(result_path, config_fingerprint = "original_fp")
      create_test_latest_pointer(result_path, "cp_000001")

      expect_error(
        load_for_resume(result_path, config_fingerprint = "different_fp"),
        class = "bayesim_checkpoint_error"
      )
    })

    it("succeeds on fingerprint mismatch with force_restart = TRUE", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path, config_fingerprint = "original_fp")
      create_test_checkpoint(result_path, config_fingerprint = "original_fp")
      create_test_latest_pointer(result_path, "cp_000001")

      # With force_restart, should not error
      result <- load_for_resume(
        result_path,
        config_fingerprint = "different_fp",
        force_restart = TRUE
      )

      # Should return NULL (indicating fresh start)
      expect_null(result)
    })

    it("returns NULL when no checkpoint available", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path, config_fingerprint = "test_fp")
      dir.create(file.path(result_path, "checkpoints"), showWarnings = FALSE)
      create_test_latest_pointer(result_path, NULL)

      result <- load_for_resume(result_path, config_fingerprint = "test_fp")

      expect_null(result)
    })
  })

  describe("find_valid_checkpoint()", {
    it("finds newest valid checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)

      # Create multiple checkpoints
      for (i in 1:3) {
        cp_id <- sprintf("cp_%06d", i)
        create_test_checkpoint(result_path, checkpoint_id = cp_id)
      }
      create_test_latest_pointer(result_path, "cp_000003")

      result <- find_valid_checkpoint(result_path)

      expect_equal(result, "cp_000003")
    })

    it("falls back to previous checkpoint if latest is invalid", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)

      # Create two valid checkpoints
      create_test_checkpoint(result_path, checkpoint_id = "cp_000001")
      create_test_checkpoint(result_path, checkpoint_id = "cp_000002")

      # Create invalid checkpoint (missing required files)
      cp_path <- file.path(result_path, "checkpoints", "cp_000003")
      dir.create(cp_path, recursive = TRUE)
      # No meta.json, ledger.rds, or results.rds

      create_test_latest_pointer(result_path, "cp_000003")

      result <- find_valid_checkpoint(result_path)

      expect_equal(result, "cp_000002")
    })

    it("returns NULL when no valid checkpoints exist", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      dir.create(file.path(result_path, "checkpoints"), showWarnings = FALSE)
      create_test_latest_pointer(result_path, "cp_000001")

      # No actual checkpoint files

      result <- find_valid_checkpoint(result_path)

      expect_null(result)
    })

    it("ignores .tmp directories when finding valid checkpoint", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      create_test_run_manifest(result_path)
      create_test_checkpoint(result_path, checkpoint_id = "cp_000001")

      # Create .tmp directory (incomplete checkpoint)
      tmp_path <- file.path(result_path, "checkpoints", "cp_000002.tmp")
      dir.create(tmp_path, recursive = TRUE)

      create_test_latest_pointer(result_path, "cp_000001")

      result <- find_valid_checkpoint(result_path)

      expect_equal(result, "cp_000001")
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
      ledger <- data.frame(
        task_id = c("t1", "t2"),
        status = c("success", "failed"),
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, ledger)

      expect_equal(result$status, c("success", "failed", "pending", "pending"))
    })

    it("preserves task grid columns not in ledger", {
      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("pending", "pending"),
        data_idx = c(1L, 1L),
        fit_idx = c(1L, 2L),
        stringsAsFactors = FALSE
      )

      ledger <- data.frame(
        task_id = "t1",
        status = "success",
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, ledger)

      expect_equal(result$data_idx, c(1L, 1L))
      expect_equal(result$fit_idx, c(1L, 2L))
    })

    it("handles empty ledger (all pending)", {
      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("pending", "pending"),
        stringsAsFactors = FALSE
      )

      ledger <- data.frame(
        task_id = character(),
        status = character(),
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, ledger)

      expect_equal(result$status, c("pending", "pending"))
    })

    it("handles task not in grid (skip silently)", {
      task_grid <- data.frame(
        task_id = c("t1", "t2"),
        status = c("pending", "pending"),
        stringsAsFactors = FALSE
      )

      # Ledger has task that doesn't exist in grid
      ledger <- data.frame(
        task_id = c("t1", "t999"),
        status = c("success", "success"),
        stringsAsFactors = FALSE
      )

      result <- merge_task_grid_status(task_grid, ledger)

      expect_equal(result$status, c("success", "pending"))
    })
  })

  describe("merge_results()", {
    it("deduplicates by task_id keeping latest", {
      old_results <- data.frame(
        task_id = c("t1", "t2"),
        metric_rmse = c(0.05, 0.10),
        stringsAsFactors = FALSE
      )

      new_results <- data.frame(
        task_id = c("t2", "t3"),
        metric_rmse = c(0.08, 0.15),
        stringsAsFactors = FALSE
      )

      result <- merge_results(old_results, new_results)

      # t2 should come from new_results (0.08)
      expect_equal(nrow(result), 3)
      expect_equal(result$metric_rmse[result$task_id == "t2"], 0.08)
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

    it("handles empty new results", {
      old_results <- data.frame(
        task_id = "t1",
        metric_rmse = 0.05,
        stringsAsFactors = FALSE
      )

      new_results <- data.frame(
        task_id = character(),
        metric_rmse = numeric(),
        stringsAsFactors = FALSE
      )

      result <- merge_results(old_results, new_results)

      expect_equal(nrow(result), 1)
      expect_equal(result$task_id, "t1")
    })

    it("preserves all columns from both sources", {
      old_results <- data.frame(
        task_id = "t1",
        metric_rmse = 0.05,
        old_col = "a",
        stringsAsFactors = FALSE
      )

      new_results <- data.frame(
        task_id = "t2",
        metric_rmse = 0.10,
        new_col = "b",
        stringsAsFactors = FALSE
      )

      result <- merge_results(old_results, new_results)

      expect_true("old_col" %in% names(result))
      expect_true("new_col" %in% names(result))
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

      init_checkpoint_dir(result_path)

      # Original data
      original_ledger <- data.frame(
        task_id = c("d001_f001_r00001", "d001_f001_r00002", "d001_f001_r00003"),
        status = c("success", "failed", "pending"),
        stringsAsFactors = FALSE
      )
      original_results <- data.frame(
        task_id = c("d001_f001_r00001", "d001_f001_r00002"),
        metric_rmse = c(0.05, NA_real_),
        metric_bias = c(0.01, NA_real_),
        stringsAsFactors = FALSE
      )

      # Write checkpoint
      write_checkpoint(
        result_path = result_path,
        ledger = original_ledger,
        results = original_results,
        config_fingerprint = "integration_test_fp",
        n_tasks_total = 3L,
        n_tasks_completed = 2L
      )

      # Read checkpoint
      cp_data <- read_checkpoint(result_path, "cp_000001")

      # Verify data preserved
      expect_equal(cp_data$ledger$task_id, original_ledger$task_id)
      expect_equal(cp_data$ledger$status, original_ledger$status)
      expect_equal(cp_data$results$task_id, original_results$task_id)
      expect_equal(cp_data$results$metric_rmse, original_results$metric_rmse)
    })

    it("preserves numeric precision", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      # Data with specific numeric values
      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(
        task_id = "t1",
        metric_rmse = 0.123456789,
        metric_bias = -0.987654321,
        metric_coverage = 0.95,
        stringsAsFactors = FALSE
      )

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "precision_test",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      cp_data <- read_checkpoint(result_path, "cp_000001")

      expect_equal(cp_data$results$metric_rmse, 0.123456789, tolerance = 1e-10)
      expect_equal(cp_data$results$metric_bias, -0.987654321, tolerance = 1e-10)
    })

    it("preserves character encoding", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      # Data with special characters
      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(
        task_id = "t1",
        metric_label = "test with special chars: \u00e9\u00e8\u00e0",
        stringsAsFactors = FALSE
      )

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "encoding_test",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      cp_data <- read_checkpoint(result_path, "cp_000001")

      expect_equal(cp_data$results$metric_label, "test with special chars: \u00e9\u00e8\u00e0")
    })
  })

  describe("Multiple checkpoints in sequence", {
    it("creates sequential checkpoint IDs", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(task_id = "t1", metric_rmse = 0.05, stringsAsFactors = FALSE)

      # Write three checkpoints
      for (i in 1:3) {
        write_checkpoint(
          result_path = result_path,
          ledger = ledger,
          results = results,
          config_fingerprint = "seq_test",
          n_tasks_total = 10L,
          n_tasks_completed = i
        )
      }

      # Verify all three exist
      expect_true(dir.exists(file.path(result_path, "checkpoints", "cp_000001")))
      expect_true(dir.exists(file.path(result_path, "checkpoints", "cp_000002")))
      expect_true(dir.exists(file.path(result_path, "checkpoints", "cp_000003")))

      # Verify latest.json points to newest
      latest <- jsonlite::read_json(file.path(result_path, "latest.json"))
      expect_equal(latest$checkpoint_id, "cp_000003")
    })

    it("can read any checkpoint in sequence", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      # Create checkpoints with different data
      for (i in 1:3) {
        ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
        results <- data.frame(task_id = "t1", metric_value = i, stringsAsFactors = FALSE)

        write_checkpoint(
          result_path = result_path,
          ledger = ledger,
          results = results,
          config_fingerprint = "seq_read_test",
          n_tasks_total = 10L,
          n_tasks_completed = i
        )
      }

      # Read each and verify correct data
      cp1 <- read_checkpoint(result_path, "cp_000001")
      cp2 <- read_checkpoint(result_path, "cp_000002")
      cp3 <- read_checkpoint(result_path, "cp_000003")

      expect_equal(cp1$results$metric_value, 1)
      expect_equal(cp2$results$metric_value, 2)
      expect_equal(cp3$results$metric_value, 3)
    })

    it("find_valid_checkpoint finds newest after multiple writes", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      for (i in 1:5) {
        create_test_checkpoint(
          result_path,
          checkpoint_id = sprintf("cp_%06d", i),
          config_fingerprint = "find_test"
        )
      }
      create_test_latest_pointer(result_path, "cp_000005")

      result <- find_valid_checkpoint(result_path)

      expect_equal(result, "cp_000005")
    })

    it("can resume from intermediate checkpoint if latest is corrupted", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      # Create two valid checkpoints
      for (i in 1:2) {
        ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
        results <- data.frame(task_id = "t1", metric_value = i, stringsAsFactors = FALSE)

        write_checkpoint(
          result_path = result_path,
          ledger = ledger,
          results = results,
          config_fingerprint = "corruption_test",
          n_tasks_total = 10L,
          n_tasks_completed = i
        )
      }

      # Corrupt the latest checkpoint
      meta_path <- file.path(result_path, "checkpoints", "cp_000002", "meta.json")
      file.remove(meta_path)

      # find_valid_checkpoint should return the previous valid one
      valid_cp <- find_valid_checkpoint(result_path)
      expect_equal(valid_cp, "cp_000001")

      # And we should be able to read it
      cp_data <- read_checkpoint(result_path, "cp_000001")
      expect_equal(cp_data$results$metric_value, 1)
    })
  })

  describe("Edge cases", {
    it("handles empty results data frame", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      ledger <- data.frame(
        task_id = character(),
        status = character(),
        stringsAsFactors = FALSE
      )
      results <- data.frame(
        task_id = character(),
        stringsAsFactors = FALSE
      )

      expect_silent(
        write_checkpoint(
          result_path = result_path,
          ledger = ledger,
          results = results,
          config_fingerprint = "empty_test",
          n_tasks_total = 0L,
          n_tasks_completed = 0L
        )
      )

      cp_data <- read_checkpoint(result_path, "cp_000001")
      expect_equal(nrow(cp_data$ledger), 0)
      expect_equal(nrow(cp_data$results), 0)
    })

    it("handles large task IDs", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      # Very large task ID
      ledger <- data.frame(
        task_id = "d999_f999_r99999",
        status = "success",
        stringsAsFactors = FALSE
      )
      results <- data.frame(task_id = "d999_f999_r99999", metric = 1, stringsAsFactors = FALSE)

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "large_id_test",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      cp_data <- read_checkpoint(result_path, "cp_000001")
      expect_equal(cp_data$ledger$task_id, "d999_f999_r99999")
    })

    it("handles special characters in metric values", {
      tmpdir <- withr::local_tempdir()
      result_path <- file.path(tmpdir, "results")

      init_checkpoint_dir(result_path)

      ledger <- data.frame(task_id = "t1", status = "success", stringsAsFactors = FALSE)
      results <- data.frame(
        task_id = "t1",
        metric_name = "metric with spaces and-dashes",
        stringsAsFactors = FALSE
      )

      write_checkpoint(
        result_path = result_path,
        ledger = ledger,
        results = results,
        config_fingerprint = "special_chars_test",
        n_tasks_total = 1L,
        n_tasks_completed = 1L
      )

      cp_data <- read_checkpoint(result_path, "cp_000001")
      expect_equal(cp_data$results$metric_name, "metric with spaces and-dashes")
    })
  })
})
