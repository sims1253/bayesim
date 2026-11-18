# Keys, identities and what callbacks see. Invariant 2 of R/study-store.R:
# everything a callback can observe is in the key of what it produces.

one_method <- function(s, fit = function(data, context) list(v = data)) {
  s <- with_methods(s, m = study_method(fit))
  with_measure(s, "v", function(artifacts, context) {
    data.frame(value = 0, seen = artifacts$fit$v)
  }, needs = "fit")
}

# S15 ------------------------------------------------------------------------------

test_that("a summary of values that are all missing has no rows and the usual columns", {
  s <- study("na", data.frame(condition_id = "a"), function(condition, context) 1)
  s <- with_methods(s, m = study_method(function(data, context) list()))
  s <- with_measure(s, "v", function(artifacts, context) {
    data.frame(target = "t", value = NA_real_)
  }, needs = "fit")
  run <- run_study(s, 3)
  out <- summarise_study(run)
  expect_equal(nrow(out), 0)
  expect_named(out, c("condition_id", "method_id", "measurement", "target", "n", "estimate", "mcse"))
  aimed <- summarise_study(run, target_mcse = 0.1)
  expect_named(aimed, c(names(out), "n_needed", "done"))
})

