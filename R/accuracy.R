#' Conduct accuracy analysis
#'
#' @description Compute ell-value, u-value, and h-value
#' @param result an object of class “peri”, usually, a result of a call to
#'   \code{\link{peri}}.
#' @param zeta A vector of variance reduction factor ([0, 1)).
#'   Default is between 0.05 and 0.8. A particular value of zeta means that we 
#'   want to reduce the variance to zeta * 100%, that is, the variance will be 
#'   multiplied by (1 - zeta). For example, if we want to reduce the variance by
#'   20%, we set zeta = 0.2. A larger value of zeta requires more observations. 
#' @example inst/examples/example-accuracy.R
#' @importFrom purrr map
#' @export
accuracy <- function(
  result,
  zeta = seq(from = 0.05, to = 0.8, by = 0.05)
) {
  zeta <- as.vector(zeta)
  n_ell <- result$data_list$n_ell
  n <- result$data_list$n_full - n_ell
  prop <- result$data_list$prop
  elss <- result$additional_info$var_estimates$elss
  output <- purrr::map(elss, \(x) .GetAllValues(zeta, n_ell, n, x, prop))

  class(output) <- c(class(output), "accuracy")
  output
}

#' Compute all metrics
#' @noRd 
.GetAllValues <- function(zeta, n_ell, n, elss, prop) {
  R_sq <- (1 - n_ell / elss) / (1 - prop)
  elss_new <- elss / (1 - zeta)
  ellvalue <- .ellvalue(zeta, n_ell, n, elss, prop, R_sq)
  uvalue <- .uvalue(zeta, n_ell, n, elss, prop, R_sq)
  hvalue <- .hvalue(zeta, n_ell, n, elss, prop, R_sq)

  all_values <- data.frame(
    zeta = zeta,
    elss = elss_new,
    ellvalue = ellvalue,
    uvalue = uvalue,
    hvalue = hvalue
  )
  
  list(
    result = all_values,
    n = n,
    n_ell = n_ell,
    elss = elss,
    r_sq = R_sq 
  )
}

#' Compute ell-value
#'
#' @noRd
.ellvalue <- function(zeta, n_ell, n, elss, prop, R_sq) {
  ellvalue <- zeta * n_ell / ((1 - R_sq) / (1 - (1 - prop) * R_sq) - zeta)
  ellvalue[zeta > (1 - elss / n)] <- NA
  ellvalue
}

#' Compute u-value
#'
#' @noRd
.uvalue <- function(zeta, n_l, n, elss, prop, R_sq) {
  v_prop <- (1 - (1 - (1 - prop) * R_sq) * (1 - zeta)) / R_sq
  uvalue <- n * (v_prop - (1 - prop)) / (1 - v_prop)
  uvalue[zeta > (prop * R_sq / (1 - (1 - prop) * R_sq))] <- NA
  uvalue
}

#' Compute h-value
#'
#' @noRd
.hvalue <- function(zeta, n_l, n, elss, prop, R_sq) {
  hvalue <- zeta * (1 - (1 - prop) * R_sq) / (1 - prop)
  hvalue[zeta > (1 - elss / n)] <- NA
  hvalue
}

#' Summary for accuracy analysis
#'
#' @param object An object of class "accuracy" returned by accuracy().
#' @param digits Number of significant digits to print. Default is 4.
#' @param ... Additional arguments (ignored).
#' @export
summary.accuracy <- function(object, digits = 4, nrows = 10, ...) {
  cat("\nAccuracy Analysis Summary\n")
  cat(strrep("-", 30), "\n")
  cat("Columns:\n")
  cat("  zeta: Variance reduction factor\n")
  cat("  elss: Effective labeled sample size\n")
  cat("  ellvalue: \u2113-value (labeled sample size needed for target variance reduction)\n")
  cat("  uvalue: u-value (unlabeled sample size needed for target variance reduction)\n")
  cat("  hvalue: h-value\n")
  cat(strrep("-", 30), "\n")
  for (i in seq_along(object)) {
    res <- object[[i]]$result
    var_name <- names(object)[i]
    if (is.null(var_name) || var_name == "") var_name <- paste0("Variable ", i)
    cat("\nVariable:", var_name, "\n")
    print(
      utils::head(
        format(
          res,
          digits = digits,
          nsmall = 2
        ),
        n = min(nrows, nrow(res))
      ),
      row.names = FALSE
    )
    cat("... (showing first", nrows, "rows)\n")
  }
  invisible(object)
}
