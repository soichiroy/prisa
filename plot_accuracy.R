#' Plot results of the accuracy analysis
#'
#' @export
#'
#' @param x an object of class “accuracy”, usually, a result of a call to \code{\link{accuracy}}.
#' @param ... other graphical parameters (see \code{\link[graphics]{par}}).
#'
plot.accuracy <- function (x, ...) {
	result <- x$result
  graphics::plot(x = result$zeta, 
                 y = result$ellvalue,
                 type = "l",
                 cex = 1.7,
                 lwd = 2.5,
                 xlab = "Variance reduction factor", 
                 ylab = paste0("\u2113", "-value"),
                 ...)
  graphics::abline(h = 0, lwd = 1.2)
  graphics::abline(v = 0, lwd = 1.2)

  graphics::plot(x = result$zeta, 
                 y = result$uvalue,
                 type = "l",
                 cex = 1.7,
                 lwd = 2.5,
                 xlab = "Variance reduction factor", 
                 ylab = "u-value",
                 ...)
  graphics::abline(h = 0, lwd = 1.2)
  graphics::abline(v = 0, lwd = 1.2)

  graphics::plot(x = result$zeta, 
                 y = result$hvalue,
                 type = "l",
                 cex = 1.7,
                 lwd = 2.5,
                 xlab = "Variance reduction factor", 
                 ylab = "h-value",
                 ...)
  graphics::abline(h = 0, lwd = 1.2)
  graphics::abline(v = 0, lwd = 1.2)
}
