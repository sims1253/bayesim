# test-resume.R
# Tests for resume pipeline: can_resume, load_for_resume, merge_task_grid_status, merge_results

library(bayesim)

test_data_gen <- function(data_spec, seed, task_ctx) {
  set.seed(seed)
  n <- data_spec$n %||% 100
  list(
    train = data.frame(x = rnorm(n), y = rnorm(n)),
    test = NULL,
    response = "y",
    true_params = c(beta = 1.0),
    vars_of_interest = "beta",
    meta = list()
  )
}

describe("can_resume()", {
  it("returns FALSE for NULL path", {
    expect_false(can_resume(NULL))
  })

  it("returns FALSE for nonexistent directory", {
    expect_false(can_resume("/tmp/bayesim_test_nonexistent_dir"))
  })
})

describe("merge_task_grid_status()", {
  it("returns fresh grid when checkpoint has no terminal tasks", {
    config <- simulation_config(
      data_grid = data.frame(n = 10),
      fit_grid = data.frame(model = "a"),
      data_generator = test_data_gen,
      seed = 42L
    )
    fresh <- create_task_grid(config)
    checkpoint_grid <- fresh
    checkpoint_grid$status <- "pending"

    merged <- merge_task_grid_status(fresh, checkpoint_grid)
    expect_equal(merged$status, rep("pending", nrow(merged)))
  })

  it("merges success status from checkpoint", {
    config <- simulation_config(
      data_grid = data.frame(n = c(10, 20)),
      fit_grid = data.frame(model = "a"),
      data_generator = test_data_gen,
      n_replicates = 2L,
      seed = 42L
    )
    fresh <- create_task_grid(config)
    checkpoint_grid <- fresh
    checkpoint_grid$status[1] <- "success"
    checkpoint_grid$status[2] <- "failed"

    merged <- merge_task_grid_status(fresh, checkpoint_grid)
    expect_equal(merged$status[1], "success")
    expect_equal(merged$status[2], "failed")
    expect_equal(merged$status[3], "pending")
  })

  it("preserves task IDs from fresh grid", {
    config <- simulation_config(
      data_grid = data.frame(n = 10),
      fit_grid = data.frame(model = "a"),
      data_generator = test_data_gen,
      seed = 42L
    )
    fresh <- create_task_grid(config)
    checkpoint_grid <- fresh
    checkpoint_grid$status[1] <- "success"

    merged <- merge_task_grid_status(fresh, checkpoint_grid)
    expect_equal(merged$task_id, fresh$task_id)
  })

  it("handles empty checkpoint grid", {
    config <- simulation_config(
      data_grid = data.frame(n = 10),
      fit_grid = data.frame(model = "a"),
      data_generator = test_data_gen,
      seed = 42L
    )
    fresh <- create_task_grid(config)
    checkpoint_grid <- fresh[0, ]

    merged <- merge_task_grid_status(fresh, checkpoint_grid)
    expect_equal(nrow(merged), nrow(fresh))
    expect_equal(merged$status, rep("pending", nrow(merged)))
  })
})

describe("merge_results()", {
  it("returns new_results when prior is NULL", {
    new_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    result <- merge_results(NULL, new_df)
    expect_equal(result, new_df)
  })

  it("returns new_results when prior is empty", {
    prior_df <- data.frame(task_id = character(), status = character())
    new_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    result <- merge_results(prior_df, new_df)
    expect_equal(nrow(result), 1)
  })

  it("returns prior when new is NULL", {
    prior_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    result <- merge_results(prior_df, NULL)
    expect_equal(result, prior_df)
  })

  it("returns prior when new is empty", {
    prior_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    new_df <- data.frame(task_id = character(), status = character())
    result <- merge_results(prior_df, new_df)
    expect_equal(nrow(result), 1)
  })

  it("combines non-overlapping results", {
    prior_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    new_df <- data.frame(task_id = "t2", status = "success", rmse = 0.3)
    result <- merge_results(prior_df, new_df)
    expect_equal(nrow(result), 2)
  })

  it("deduplicates overlapping task_ids", {
    prior_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    new_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    result <- merge_results(prior_df, new_df)
    expect_equal(nrow(result), 1)
  })

  it("errors on conflicting duplicate task_ids", {
    prior_df <- data.frame(task_id = "t1", status = "success", rmse = 0.5)
    new_df <- data.frame(task_id = "t1", status = "success", rmse = 0.3)
    expect_error(
      merge_results(prior_df, new_df),
      "Conflicting duplicate"
    )
  })
})

describe("get_resume_summary()", {
  it("returns NULL when cannot resume", {
    expect_null(get_resume_summary("/tmp/bayesim_test_nonexistent"))
  })
})

describe("format_resume_summary()", {
  it("returns correct format for valid summary", {
    summary <- list(checkpoint_id = 1L, n_total = 100L, n_completed = 50L, n_pending = 50L)
    result <- format_resume_summary(summary)
    expect_true(grepl("Checkpoint 1", result))
    expect_true(grepl("50/100", result))
    expect_true(grepl("50 pending", result))
  })

  it("returns descriptive message for NULL input", {
    expect_equal(format_resume_summary(NULL), "No resumable state found")
  })
})
