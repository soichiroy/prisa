#' Compute ell-value, u-value, and h-value
#' 
#' @param result an object of class “MLcovar”, usually, a result of a call to \code{\link{MLcovar}}.
#' @param zeta A vector of variance reduction factor ([0, 1]).
#'   Default is between 0.05 and 0.8.
#' @example examples/example-accuracy.R
#' @export
accuracy <- function(
  result,
  zeta = seq(from = 0.05, to = 0.8, by = 0.05)
) {

  zeta <- as.vector(zeta)
  n_ell <- result$data_list$n_ell
  n <- result$data_list$n_full - n_ell
  elss <- result$additional_info$var_estimates$elss
  prop <- result$data_list$prop
  R_sq <- (1 - n_ell / elss) / (1 - prop)
  elss_new <- elss / (1 - zeta)
  ellvalue <- .ellvalue(zeta, n_ell, n, elss, prop, R_sq)
  uvalue <- .uvalue(zeta, n_ell, n, elss, prop, R_sq)
  hvalue <- .hvalue(zeta, n_ell, n, elss, prop, R_sq)

  elluhvalues <- data.frame(
    zeta = zeta,
    elss = elss_new,
    ellvalue = ellvalue,
    uvalue = uvalue,
    hvalue = hvalue
  )
  output <- list(
    result = elluhvalues,
    n = n,
    n_ell = n_ell,
    elss = elss,
    r_sq = R_sq
  )

  class(output) <- c(class(output), "accuracy")
  output
}

#' Compute ell-value
#' 
#' @noRd 
.ellvalue <- function (zeta, n_ell, n, elss, prop, R_sq) {
  ellvalue <- zeta * n_ell / ((1 - R_sq) / (1 - (1 - prop) * R_sq) - zeta)
  ellvalue[zeta > (1 - elss / n)] <- NA
  ellvalue
}

#' Compute u-value
#' 
#' @noRd 
.uvalue <- function (zeta, n_l, n, elss, prop, R_sq) {
  v_prop <- (1 - (1 - (1 - prop) * R_sq) * (1 - zeta)) / R_sq
  uvalue <- n * (v_prop - (1 - prop)) / (1 - v_prop)
  uvalue[zeta > (prop * R_sq / (1 - (1 - prop) * R_sq))] <- NA
  uvalue
}

#' Compute h-value
#' 
#' @noRd 
.hvalue <- function (zeta, n_l, n, elss, prop, R_sq) {
  hvalue <- zeta * (1 - (1 - prop) * R_sq) / (1 - prop)
  hvalue[zeta > (1 - elss / n)] <- NA
  hvalue
}
