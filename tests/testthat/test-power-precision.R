# I5: n_replicates_for_target() — invert MCSE formulas for planning.
library(bayesim)

describe("n_replicates_for_target", {
  it("coverage default p=0.5: n = ceil(0.25 / mcse^2)", {
    # target_mcse = 0.03 -> 0.25 / 0.0009 = 277.7... -> 278
    expect_equal(n_replicates_for_target(0.03, "coverage"), 278L)
  })

  it("coverage p=0.9: n = ceil(0.09 / mcse^2)", {
    # 0.09 / 0.0009 = 100
    expect_equal(n_replicates_for_target(0.03, "coverage", p = 0.9), 100L)
  })

  it("continuous: n = ceil((sd / mcse)^2)", {
    # (0.5 / 0.05)^2 = 100
    expect_equal(
      n_replicates_for_target(0.05, "continuous", assumed_sd = 0.5),
      100L
    )
  })

  it("returns an integer scalar", {
    expect_true(is.integer(n_replicates_for_target(0.01, "coverage")))
    expect_length(n_replicates_for_target(0.01, "coverage"), 1L)
  })

  it("metric_type defaults to coverage", {
    expect_equal(
      n_replicates_for_target(0.03),
      n_replicates_for_target(0.03, "coverage")
    )
  })

  it("errors on bad target_mcse", {
    expect_error(
      n_replicates_for_target(0, "coverage"),
      class = "bayesim_config_error"
    )
    expect_error(
      n_replicates_for_target(-1, "coverage"),
      class = "bayesim_config_error"
    )
    expect_error(
      n_replicates_for_target(NA_real_, "coverage"),
      class = "bayesim_config_error"
    )
    expect_error(
      n_replicates_for_target(c(0.01, 0.02), "coverage"),
      class = "bayesim_config_error"
    )
  })

  it("errors on bad p for coverage", {
    expect_error(
      n_replicates_for_target(0.01, "coverage", p = 1.5),
      class = "bayesim_config_error"
    )
    expect_error(
      n_replicates_for_target(0.01, "coverage", p = -0.1),
      class = "bayesim_config_error"
    )
  })

  it("errors when assumed_sd missing/invalid for continuous", {
    expect_error(
      n_replicates_for_target(0.01, "continuous"),
      class = "bayesim_config_error"
    )
    expect_error(
      n_replicates_for_target(0.01, "continuous", assumed_sd = 0),
      class = "bayesim_config_error"
    )
    expect_error(
      n_replicates_for_target(0.01, "continuous", assumed_sd = -1),
      class = "bayesim_config_error"
    )
  })

  it("errors on unknown metric_type", {
    expect_error(n_replicates_for_target(0.01, "nope"), class = "simpleError")
  })
})
