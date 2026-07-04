# Workstream I4: report() renders a Quarto HTML report.

skip_on_cran()

.gen <- function(data_spec, task_ctx) {
  n <- data_spec$n %||% 50L
  beta <- data_spec$beta %||% 0.5
  x <- stats::rnorm(n)
  y <- 2 + beta * x + stats::rnorm(n)
  list(
    train = data.frame(y = y, x = x),
    test = NULL,
    response = "y",
    true_params = c(Intercept = 2, x = beta, sigma = 1),
    vars_of_interest = c("Intercept", "x", "sigma"),
    meta = list()
  )
}

describe("report()", {
  it("renders an HTML report for a small LinearRegressionFitter study", {
    skip_if_not(
      requireNamespace("quarto", quietly = TRUE) &&
        (nzchar(Sys.which("quarto")) || !is.null(quarto::quarto_path())),
      "quarto CLI not available"
    )

    config <- simulation_config(
      data_grid = data.frame(n = 50L, beta = 0.5),
      fit_grid = data.frame(model = "lm"),
      data_generator = .gen,
      fitter = LinearRegressionFitter(n_draws = 100L),
      metrics = list(posterior_summary_metric()),
      n_replicates = 3L,
      seed = 7L
    )
    result <- run_simulation(config, resume = "never", progress = FALSE)

    out <- tempfile(fileext = ".html")
    res <- report(result, output_file = out, open = FALSE)
    expect_true(file.exists(res))
    expect_gt(file.info(res)$size, 0)
  })
})
