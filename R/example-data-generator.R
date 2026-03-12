#' @keywords internal
bayesim_example_data_generator <- function(data_spec, seed, task_ctx) {
  withr::with_seed(seed, {
    n <- as.integer(data_spec$n %||% 10L)
    beta <- as.numeric(data_spec$beta %||% 1)
    sigma <- as.numeric(data_spec$sigma %||% 1)

    x <- stats::rnorm(n)
    y <- beta * x + stats::rnorm(n, sd = sigma)

    list(
      train = data.frame(y = y, x = x),
      test = NULL,
      response = "y",
      true_params = c(beta = beta, sigma = sigma),
      vars_of_interest = c("beta", "sigma"),
      references = c(beta = 0, sigma = 1),
      meta = list(task_id = task_ctx$task_id %||% NA_character_)
    )
  })
}
