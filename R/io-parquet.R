# =============================================================================
# Optional parquet summary output (Workstream I8)
# =============================================================================
#
# The wide rds summary becomes a bottleneck for very large studies. This file
# adds an OPTIONAL parquet sidecar for the summary, written alongside the
# canonical rds checkpoint. The rds checkpoint remains the resume artifact;
# the parquet file is for downstream (pandas/arrow/polars) consumption.
#
# `nanoparquet` is a lightweight parquet reader/writer with no system
# dependencies, declared in Suggests.

#' Write a simulation summary to parquet
#'
#' Writes `results_df` (a data frame/tibble summary produced by
#' [build_simulation_result()]) to `path` as a parquet file using
#' `nanoparquet`. The data frame is coerced to a plain data.frame first so
#' that tibble/vctrs columns parquet cannot represent are flattened.
#'
#' Requires the suggested `nanoparquet` package; if it is not installed,
#' [rlang::check_installed()] prompts (or errors) accordingly.
#'
#' @param results_df A data.frame or tibble summary to write.
#' @param path Character scalar. Destination parquet file path.
#' @return Invisible `path` (invisibly), the path written to.
#' @keywords internal
write_results_parquet <- function(results_df, path) {
  rlang::check_installed("nanoparquet")
  nanoparquet::write_parquet(as.data.frame(results_df), path)
  invisible(path)
}

#' Read a simulation summary file
#'
#' Reads a simulation summary written by [run_simulation()] (or
#' [write_results_parquet()]). Dispatches on file extension:
#' `.parquet` files are read with `nanoparquet::read_parquet` (requires the
#' suggested `nanoparquet` package); `.rds` files are read with
#' [base::readRDS()].
#'
#' @param path Character scalar. Path to a `.parquet` or `.rds` summary file.
#' @return A data frame of the summary.
#' @export
#' @examples
#' \dontrun{
#' summary_df <- read_summary("results/summary.parquet")
#' summary_df <- read_summary("results/summary.rds")
#' }
read_summary <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    cli::cli_abort("{.arg path} must be a single character string")
  }
  if (!file.exists(path)) {
    cli::cli_abort("File does not exist: {.file {path}}")
  }

  ext <- tolower(tools::file_ext(path))
  if (identical(ext, "parquet")) {
    rlang::check_installed("nanoparquet")
    nanoparquet::read_parquet(path)
  } else if (identical(ext, "rds")) {
    readRDS(path)
  } else {
    cli::cli_abort(c(
      "Unsupported summary file extension: {.val {ext}}",
      i = "{.fn read_summary} supports {.val parquet} and {.val rds}."
    ))
  }
}
