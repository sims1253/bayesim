#' Title
#'
#' @param dist
#' @param bounds
#' @param pars
#' @param prefix
#' @param parnames
#' @param package
#' @param user_theme
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_dist <- function(dist, bounds, pars, prefix = "d", parnames = NULL,
                      package = NULL, user_theme = ggthemes::theme_tufte, ...) {
  `%>%` <- dplyr::`%>%`
  pos <- -1
  if (!is.null(package)) {
    pos <- asNamespace(package)
  }
  ddist <- get(paste0(prefix, dist), pos = pos, mode = "function")
  df <- data.frame(x = seq(bounds[1], bounds[2], 0.001))
  if (!is.null(parnames)) {
    parnames <- paste0(parnames, " = ")
  }
  cnames <- rep(NA, length(pars))
  for (i in seq_along(pars)) {
    tmp <- do.call(ddist, c(list(df$x), pars[[i]], list(...)))
    cnames[i] <- paste0("$", parnames, pars[[i]], "$", collapse = ", ")
    df[paste0(parnames, pars[[i]], collapse = ", ")] <- tmp
  }
  df <- df %>%
    tidyr::gather("pars", "dens", -x) %>%
    dplyr::mutate(pars = factor(pars, unique(pars)))
  gg <- ggplot2::ggplot(df, ggplot2::aes(x, dens, color = pars)) +
    user_theme() +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_color_viridis_d(labels = unname(latex2exp::TeX(cnames))) +
    ggplot2::labs(x = "x", y = "", color = "") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10)
    )
  if (prefix == "p") {
    gg <- gg +
      ggplot2::scale_y_continuous(breaks = c(0, 0.5, 1)) +
      ggplot2::theme(axis.ticks.y = ggplot2::element_line(),
            axis.text.y = ggplot2::element_text(),
            axis.line.y = ggplot2::element_line())
  } else if (prefix == "q") {
    gg <- gg +
      ggplot2::scale_y_continuous() +
      ggplot2::theme(axis.ticks.y = ggplot2::element_line(),
            axis.text.y = ggplot2::element_text(),
            axis.line.y = ggplot2::element_line())
  }
  return(gg)
}
