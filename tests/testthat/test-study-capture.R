# Capture: what a callback takes to its worker, and what its identity
# counts. One function, .capture() in R/study-capture.R, decides both, so a
# binding is in both or in neither.

identity_of <- function(f) bayesim:::.capture_function(f, "f")$id
transport <- function(f) bayesim:::.capture_function(f, "f")$fn

# A callback that returns `x`.
returning <- function(x) {
  force(x)
  function() x
}

one_method <- function(fit) {
  s <- study("capture", data.frame(condition_id = "a"), function(condition, context) 1)
  with_methods(s, m = study_method(fit))
}

# Bind values in the global environment for one test, as a user's script
# does. A worker's global environment does not hold them.
local_globals <- function(..., env = parent.frame()) {
  values <- list(...)
  list2env(values, globalenv())
  withr::defer(rm(list = names(values), envir = globalenv()), envir = env)
}

in_global <- function(expr) eval(substitute(expr), globalenv())

# Call a transported callback in a new R process.
on_worker <- function(fn, ...) {
  callr::r(function(fn, args) do.call(fn, args), list(fn, list(...)))
}

# S20 --------------------------------------------------------------------------------

test_that("an S4 object counts with its data part", {
  make <- function(v) returning(asS4(structure(v, tag = "a")))
  expect_false(identical(identity_of(make(1:3)), identity_of(make(4:6))))
})

test_that("the order of attributes is part of the identity", {
  ab <- structure(1, a = 1, b = 2)
  ba <- structure(1, b = 2, a = 1)
  expect_false(identical(identity_of(returning(ab)), identity_of(returning(ba))))
})

test_that("planning refuses an active binding and names it and the callback", {
  e <- new.env()
  makeActiveBinding("now", function() Sys.time(), e)
  fit <- local(function(data, context) list(t = now), envir = e)
  expect_error(
    plan_study(one_method(fit), 1),
    "Method \"m\" \\(fit\\) captures .now.*active binding"
  )
})

test_that("planning refuses a promise that is not forced and names it and the callback", {
  make <- function(k) function(data, context) list(k = k)
  # The callback would read `k` when it first runs, after `k` has changed.
  k <- 1
  fit <- make(k)
  k <- 2
  expect_error(
    plan_study(one_method(fit), 1),
    "Method \"m\" \\(fit\\) captures .k.*not forced"
  )
  # A promise of a value gives that value whenever it is forced.
  expect_s3_class(plan_study(one_method(make(1)), 1), "bayesim_study_plan")
  # A forced argument travels.
  make_forced <- function(k) {
    force(k)
    function(data, context) list(k = k)
  }
  expect_s3_class(plan_study(one_method(make_forced(1)), 1), "bayesim_study_plan")
})

# S21 --------------------------------------------------------------------------------

test_that("a formula whose environment holds the formula ends the walk", {
  e <- new.env(parent = baseenv())
  f <- as.formula("~ f", env = e)
  e$f <- f
  cb <- function() f
  expect_type(identity_of(cb), "character")
  # The copy keeps the cycle.
  fo <- get("f", environment(transport(cb)))
  expect_identical(environment(fo)$f, fo)
})

# S22 --------------------------------------------------------------------------------

test_that("a closure in a list travels with what it reads", {
  skip_if_not_installed("callr")
  local_globals(x = 5)
  local_globals(helpers = in_global(list(f = function() x)))
  fit <- in_global(function(data, context) list(v = helpers$f()))
  expect_equal(on_worker(transport(fit), NULL, NULL)$v, 5)
  # What travels is what counts.
  before <- identity_of(fit)
  assign("x", 6, globalenv())
  expect_false(identical(identity_of(fit), before))
})

# S23 --------------------------------------------------------------------------------

test_that("keys and column orders do not depend on the collation locale", {
  skip_if_not_installed("callr")
  in_locale <- function(locale) {
    callr::r(function() {
      B <- 1
      a <- 2
      `_z` <- 3
      Z <- 4
      reads <- function() c(B, a, `_z`, Z)
      s <- bayesim::study(
        "locale",
        data.frame(condition_id = "c", B = 1, a = 2),
        function(condition, context) reads()
      )
      s <- bayesim::with_methods(s, m = bayesim::study_method(
        function(data, context) list(v = data),
        settings = list(B = 1, a = 2)
      ))
      s <- bayesim::with_measure(s, "B", function(artifacts, context) {
        data.frame(value = 1)
      }, needs = "fit")
      s <- bayesim::with_measure(s, "a", function(artifacts, context) {
        data.frame(value = 2)
      }, needs = "fit")
      setup <- bayesim:::.sg_setup(s, 1, 1)
      run <- bayesim::run_study(s, 1)
      list(
        keys = c(setup$datasets$dataset_id, setup$fits$record),
        columns = names(run$measurements)
      )
    }, env = c(callr::rcmd_safe_env(), LC_ALL = locale))
  }
  expect_identical(in_locale("C.UTF-8"), in_locale("C"))
})

# S31 --------------------------------------------------------------------------------

test_that("a replacement of a free variable reads and counts the variable", {
  make <- function(b) {
    x <- list(b = b)
    function() {
      x$a <- 1
      x$b
    }
  }
  expect_false(identical(identity_of(make(2)), identity_of(make(3))))
  skip_if_not_installed("callr")
  local_globals(x = list(b = 2))
  cb <- in_global(function() {
    x$a <- 1
    x$b
  })
  expect_equal(on_worker(transport(cb)), 2)
})

test_that("planning refuses the dots of an enclosing function", {
  make <- function(...) function(data, context) list(...)
  expect_error(
    plan_study(one_method(make(1)), 1),
    "Method \"m\" \\(fit\\) captures .*\\.\\.\\..*enclosing function"
  )
})
