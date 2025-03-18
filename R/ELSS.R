#' Compute Effective Label Sample Size (ELSS)
#' @noRd 
.ComputeELSS <- function(n_ell, coef_estimates) {
  var_unadjusted <- diag(coef_estimates$vcov)[1]
  n_elss <- n_ell * var_unadjusted / coef_estimates$var
  return(n_elss)
}
