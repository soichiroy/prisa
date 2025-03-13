

#' Compute Effective Label Sample Size (ELSS)
#' @noRd 
.ComputeELSS <- function(n_ell, coef_estimates) {
  var_unadjusted <- diag(coef_estimates$vcov)[1]
  delta <- 1 - coef_estimates$var / var_unadjusted
  n_elss <- n_ell / (1 - delta) 
  return(n_elss)
}