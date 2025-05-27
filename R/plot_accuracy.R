#' Plot results of the accuracy analysis
#'
#' @export
#' @import ggplot2
#' @import patchwork
#' @importFrom tidyr drop_na
#'
#' @param x an object of class “accuracy”, usually, a result of a call to \code{\link{accuracy}}.
#' @param ... other graphical parameters (see \code{\link[graphics]{par}}).
#'
plot.accuracy <- function(x, ...) {
  p_list <- purrr::map(x, function(i) {
    res <- i$result
    plot_accuracy_single_param(res)
  })
  p_list
}

plot_accuracy_single_param <- function(df) {
  zeta_min <- min(df$zeta, na.rm = TRUE)
  zeta_max <- max(df$zeta, na.rm = TRUE)
  plot_theme <- theme_bw() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  p1 <- df %>%
    drop_na(.data$zeta, .data$ellvalue) %>%
    ggplot(aes(x = zeta, y = ellvalue)) +
    geom_line(linewidth = 1.1) +
    scale_x_continuous(
      limits = c(zeta_min, zeta_max),
      labels = scales::percent_format(),
    ) +
    labs(x = "Variance reduction factor", y = paste0("\u2113", "-value")) +
    plot_theme
  p2 <- df %>%
    drop_na(.data$zeta, .data$uvalue) %>%
    ggplot(aes(x = zeta, y = uvalue)) +
    geom_line(linewidth = 1.1) +
    scale_x_continuous(
      limits = c(zeta_min, zeta_max),
      labels = scales::percent_format(),
    ) +
    labs(x = "Variance reduction factor", y = "u-value") +
    plot_theme

  p3 <- df %>%
    drop_na(.data$zeta, .data$hvalue) %>%
    ggplot(aes(x = zeta, y = hvalue)) +
    geom_line(linewidth = 1.1) +
    scale_x_continuous(
      limits = c(zeta_min, zeta_max),
      labels = scales::percent_format(),
    ) +
    labs(x = "Variance reduction factor", y = "h-value") +
    plot_theme

  # Display all three plots together
  p <- (p1 / p2 / p3)
  return(p)

}