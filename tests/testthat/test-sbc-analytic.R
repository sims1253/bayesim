# Golden workflow 6: analytic SBC acceptance.
#
# Exact conjugate prior-predictive simulation-based calibration (Talts et al.
# 2018) using LinearRegressionFitter: the generator draws (Intercept, x,
# sigma) from the same Normal-Inverse-Gamma prior the fitter updates, so the
# posterior is the exact conjugate updater and SBC must pass by construction.
# Executable port of vignettes/sbc-and-calibration.Rmd with a modest budget
# (100 replicates, 400 draws) that keeps the check inside the fast PR gate.
library(bayesim)

.sbc_gen <- function(data_spec, task_ctx) {
  n <- data_spec$n
  # Exact NIG prior used by the fitter below:
  #   sigma^2 ~ Inv-Gamma(2, 1); beta | sigma^2 ~ N(0, sigma^2 / 0.25)
  sigma <- sqrt(1 / stats::rgamma(1, shape = 2, rate = 1))
  intercept <- stats::rnorm(1, mean = 0, sd = sigma / sqrt(0.25))
  slope <- stats::rnorm(1, mean = 0, sd = sigma / sqrt(0.25))
  x <- stats::rnorm(n)
  y <- intercept + slope * x + stats::rnorm(n, sd = sigma)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(Intercept = intercept, x = slope, sigma = sigma),
    vars_of_interest = c("Intercept", "x", "sigma")
  )
}

.sbc_run <- function(seed = 7L) {
  config <- simulation_config(
    data_grid = data.frame(n = 30L),
    fit_grid = data.frame(model = "lm"),
    data_generator = .sbc_gen,
    fitter = LinearRegressionFitter(
      n_draws = 400L,
      prior_mean = 0,
      prior_precision = 0.25,
      a0 = 2,
      b0 = 1
    ),
    metrics = list(
      rank_metric(thin = "auto"),
      posterior_summary_metric()
    ),
    n_replicates = 100L,
    seed = seed
  )
  run_simulation(config, progress = FALSE)
}

describe("golden workflow: analytic SBC acceptance", {
  it("records one bounded integer rank per task for every parameter", {
    result <- .sbc_run()
    summ <- result$summary

    expect_equal(sum(summ$status == "success"), 100L)
    for (param in c("Intercept", "x", "sigma")) {
      col <- paste0("rank__by_param__", param)
      expect_true(col %in% names(summ), info = col)
      ranks <- summ[[col]]
      expect_true(all(ranks == floor(ranks)))
      expect_true(all(ranks >= 0 & ranks <= 400L), info = param)
      expect_false(anyNA(ranks))
    }
    # Retained rank counts: thin="auto" keeps stride 1 for these i.i.d.
    # conjugate draws, so S = 400 and n_ranks = S + 1.
    expect_true(all(summ$rank__n_draws == 400L))
    expect_true(all(summ$rank__n_ranks == 401L))
  })

  it("passes the uniformity check for the conjugate exact posterior", {
    ranks <- sbc_ranks(.sbc_run())

    expect_equal(nrow(ranks), 300L) # 3 params x 100 replicates
    expect_setequal(ranks$param, c("Intercept", "x", "sigma"))

    S <- max(ranks$n_ranks) - 1L
    expect_equal(S, 400L)
    # Ranks are integers in [0, S] for every parameter.
    expect_true(all(ranks$rank >= 0 & ranks$rank <= S))
    expect_true(all(ranks$rank == floor(ranks$rank)))

    # Under uniformity E[rank] = S/2; the standard error of the mean rank is
    # sqrt((S^2 - 1) / 12 / L). A tolerance of ~5 SE is generous but stable
    # for an exact conjugate updater at L = 100 replicates.
    for (param in c("Intercept", "x", "sigma")) {
      rr <- ranks$rank[ranks$param == param]
      se <- sqrt(((S + 1)^2 - 1) / 12 / length(rr))
      expect_lt(abs(mean(rr) - S / 2), 5 * se)
    }
  })

  it("is deterministic across runs with the same seed", {
    run1 <- .sbc_run(seed = 7L)
    run2 <- .sbc_run(seed = 7L)

    # Compare the scientific output; wall-clock timing columns legitimately
    # differ between runs.
    rank_cols <- grep("^rank__", names(run1$summary), value = TRUE)
    est_cols <- grep(
      "^posterior_summary__|^truth__",
      names(run1$summary),
      value = TRUE
    )
    for (col in c(rank_cols, est_cols)) {
      expect_equal(run1$summary[[col]], run2$summary[[col]], info = col)
    }
    expect_equal(run1$summary$task_id, run2$summary$task_id)
  })

  it("feeds the rank ECDF machinery (ggplot2)", {
    skip_if_not_installed("ggplot2")
    ranks <- sbc_ranks(.sbc_run())

    expect_s3_class(plot_rank_ecdf(ranks, alpha = 0.95), "ggplot")
    expect_s3_class(plot_rank_ecdf(ranks, alpha = 0.99), "ggplot")
    expect_s3_class(plot_rank_hist(ranks), "ggplot")
  })
})
