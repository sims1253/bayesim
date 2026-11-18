# Reads the source tree: the installed Rd database drops \docType, so this
# runs only where the package sources sit next to the tests.
test_that("exactly one topic documents the package", {
  root <- testthat::test_path("..", "..")
  skip_if_not(
    file.exists(file.path(root, "DESCRIPTION")) &&
      dir.exists(file.path(root, "R")) &&
      dir.exists(file.path(root, "man")),
    "package sources not available"
  )

  r_files <- list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)
  declares <- vapply(r_files, function(f) {
    any(trimws(readLines(f, warn = FALSE)) == "\"_PACKAGE\"")
  }, logical(1))
  expect_identical(basename(r_files[declares]), "bayesim-package.R")

  rd_files <- list.files(file.path(root, "man"), "\\.Rd$", full.names = TRUE)
  package_rd <- vapply(rd_files, function(f) {
    any(grepl("\\docType{package}", readLines(f, warn = FALSE), fixed = TRUE))
  }, logical(1))
  expect_identical(basename(rd_files[package_rd]), "bayesim-package.Rd")
})
