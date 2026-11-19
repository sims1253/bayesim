# Study summaries ---------------------------------------------------------------------

#' Summarise study results with Monte Carlo standard errors
#'
#' `summarise_study()` averages measurement or comparison values over
#' replicates, by condition, method, measurement and target. Each estimate
#' comes with its Monte Carlo standard error (MCSE).
#'
#' For a measurement whose summary is `"rmse"` (see [measure_rmse()]), the values are
#' squared errors; the estimate is the root of their mean, and its MCSE comes
#' from the delta method. For other measurements, the estimate is the mean.
#'
#' Comparisons are summarised per condition. When a comparison returns one
#' paired difference per dataset, as [compare_paired()] does, the MCSE is
#' the paired MCSE: the standard deviation of the differences divided by the
#' square root of their number.
#'
#' @param x A [run_study()] or [assess_study()] result.
#' @param what `"measurements"` or `"comparisons"`.
#' @param by Further columns to group by.
#' @param target_mcse Optional MCSE to aim for. The summary then adds
#'   `n_needed`, an estimate of the number of replicates that gives this
#'   MCSE, and `done`, which is `TRUE` when the MCSE is already at most
#'   `target_mcse`.
#'
#' @return A data frame with the group columns and `n`, `estimate` and
#'   `mcse`.
#' @export
summarise_study <- function(
  x,
  what = c("measurements", "comparisons"),
  by = NULL,
  target_mcse = NULL
) {
  if (!inherits(x, c("bayesim_study_run", "bayesim_study_assessment"))) {
    cli::cli_abort(
      "{.arg x} must be a {.fn run_study} or {.fn assess_study} result."
    )
  }
  what <- match.arg(what)
  if (!is.null(target_mcse) &&
      (!is.numeric(target_mcse) || length(target_mcse) != 1L || !(target_mcse > 0))) {
    cli::cli_abort("{.arg target_mcse} must be one positive number.")
  }
  rows <- x[[what]]
  if (!nrow(rows)) {
    return(data.frame())
  }
  if (!"value" %in% names(rows) || !is.numeric(rows$value)) {
    cli::cli_abort("The {what} need a numeric {.field value} column.")
  }
  if (what == "measurements") {
    keys <- c("condition_id", "method_id", "measurement")
    kinds <- vapply(x$study$measures, function(m) m$summary %||% "mean", "")
    rows$.summary <- unname(kinds[rows$measurement])
  } else {
    if (!"condition_id" %in% names(rows) && "dataset_id" %in% names(rows)) {
      lookup <- unique(x$attempts[c("dataset_id", "condition_id")])
      rows$condition_id <- lookup$condition_id[match(rows$dataset_id, lookup$dataset_id)]
    }
    keys <- intersect(c("condition_id", "comparison", "method_id", "reference"), names(rows))
    rows$.summary <- "mean"
  }
  keys <- unique(c(keys, intersect("target", names(rows)), by))
  unknown <- setdiff(keys, names(rows))
  if (length(unknown)) {
    cli::cli_abort("Unknown {.arg by} columns: {paste(unknown, collapse = ', ')}.")
  }
  rows <- rows[!is.na(rows$value), , drop = FALSE]
  g <- group_ids(rows, keys)
  first <- match(unique(g), g)
  out <- rows[first, keys, drop = FALSE]
  stats <- lapply(split(seq_len(nrow(rows)), factor(g, levels = unique(g))), function(i) {
    v <- rows$value[i]
    n <- length(v)
    m <- mean(v)
    se <- if (n > 1L) stats::sd(v) / sqrt(n) else NA_real_
    if (identical(rows$.summary[[i[[1L]]]], "rmse")) {
      est <- sqrt(m)
      c(n, est, if (est > 0) se / (2 * est) else NA_real_)
    } else {
      c(n, m, se)
    }
  })
  stats <- do.call(rbind, c(list(matrix(numeric(), 0L, 3L)), stats))
  out$n <- as.integer(stats[, 1L])
  out$estimate <- stats[, 2L]
  out$mcse <- stats[, 3L]
  if (!is.null(target_mcse)) {
    out$n_needed <- as.integer(ceiling(out$n * (out$mcse / target_mcse)^2))
    out$done <- !is.na(out$mcse) & out$mcse <= target_mcse
  }
  rownames(out) <- NULL
  out
}
