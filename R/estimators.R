#' @title Estimators for the main and proxy models
#'
#' @noRd 
#' @inheritParams .PredictionErrorRobustInference 
#' @param use_label_only Boolean. If TRUE, only the labeled set estimates are
#'   returned. Default is FALSE. 
.GetPointEstimates <- function(
    main_model,
    proxy_model,
    data_list,
    args_main_model,
    args_proxy_model,
    use_label_only = FALSE) {

  # Unbiased estimator
  tau_ell <-
    do.call(main_model, c(list(data_list$dat_labeled), args_main_model))

  # Biased estimators
  delta_ell <- do.call(
    proxy_model, c(list(data_list$dat_labeled), args_proxy_model)
  )

  n_main_estimates <- length(tau_ell)
  n_proxy_estimates <- length(delta_ell)
  # Return estimates based on the labeled set
  if (isTRUE(use_label_only)) {
    # Return estimates
    estimates <- list(
      tau_ell = tau_ell,
      delta_ell = delta_ell,
      n_main_estimates = n_main_estimates,
      n_proxy_estimates = n_proxy_estimates
    )
    return(estimates)
  }

  # Proxy estimators on the full data
  delta_full <- do.call(
    proxy_model, c(list(data_list$dat_full), args_proxy_model)
  )

  # Create control variate estimators
  delta_diff <- as.vector(delta_ell) - as.vector(delta_full)
  
  # Return estimates 
  estimates <- list(
    tau_ell = tau_ell,
    delta_diff = delta_diff,
    n_main_estimates = n_main_estimates,
    n_proxy_estimates = n_proxy_estimates
  )

  return(estimates)
}

#' Estimate optimal coefficients
#' 
#' @param cov_estimates A list of bootstrap variance-covariance estimates.
#' @noRd
.EstimateOptimalCoefficients <- function(cov_estimates) {
  # Variance covariance matrix of the proxy estimator
  vcov_full <- cov_estimates$vcov_full

  if (!is.null(vcov_full)) {
    # When the vcov_full is available, estimate the coefficient using the 
    # full formula: A* = Cov(tau_main_ell, tau_proxy_diff) / V(tau_proxy_diff)
    coef_estimates <- solve(vcov_full, cov_estimates$vcov_main_diff)
    return(coef_estimates)
  }

  #
  # When use_full is FALSE, use the labeled set estimator to estimate the
  # coefficients.
  #

  n_main_estimates <- cov_estimates$n_main_estimates
  vcov_labeled <- cov_estimates$vcov_labeled
  # Covariance between the main and proxy model estimates
  #  * [1,1] element is the variance of the main unbiased estimator 
  idx_main <- 1:n_main_estimates
  cov_main_proxy <- vcov_labeled[-idx_main, idx_main, drop = FALSE]

  # When the vcov_full is not available, use the labeled set estimator to
  # estimate the coefficients.
  #  A* = COV(tau_proxy_ell, tau_main_ell) / V(tau_proxy_ell)
  coef_estimates <-
    solve(vcov_labeled[-idx_main, -idx_main, drop = FALSE], cov_main_proxy)
  coef_estimates
}

.CombineEstimates <- function(point_estimate, coef_estimates) {
  tau_ell <- point_estimate$tau_ell
  delta_diff <- point_estimate$delta_diff

  # Proposed estimator: unbiased + coef * (cv_estimators)
  combined_estimate <- tau_ell - as.vector(t(coef_estimates) %*% delta_diff) 
  combined_estimate
}

#' Estimate the variance of the proposed estimator
#' 
#' @noRd
.EstimateVariance <- function(
    cov_estimates,
    coef_estimates,
    prop,
    n_ell,
    options) {

  n_main_estimates <- cov_estimates$n_main_estimates

  # Variance covariance of the main and proxy estimators
  vcov_labeled <- cov_estimates$vcov_labeled

  # Variance of the original estimator
  var_tau_ell <- diag(vcov_labeled)[1:n_main_estimates]

  # vcov_main_diff is NULL when use_full is FALSE.
  if (isFALSE(options$use_full)) {
    n_main_estimates <- cov_estimates$n_main_estimates
    idx_main <- 1:n_main_estimates
    cov_main_proxy <- vcov_labeled[idx_main, -idx_main, drop = FALSE]
    # Variance of the proposed estimator:
    #   Additional (1 - prop) scaling because cov_main_proxy is based on the
    #   labeled data.
    var_est <- var_tau_ell -
      (1 - prop) * diag(cov_main_proxy %*% coef_estimates) 
  } else {
    # Use the full formula to estimate the variance
    var_est <- var_tau_ell - 
      diag(t(cov_estimates$vcov_main_diff) %*% coef_estimates)
  }

  # Check the variance estimate
  var_theoretical_limit <- prop * var_tau_ell
  if (any(var_est < var_theoretical_limit) && !options$debug_mode) {
    stop(
      "Failed variance estimation. Please use use_full = TRUE option instead."
    )
  }

  if (any(var_est < var_theoretical_limit) && options$debug_mode) {
    warning(
      "Estimated variance is too small.", 
      "Using the theoretical limit instead.", 
      "Only for debugging."
    )
    var_est <- var_theoretical_limit
  }
  

  # Compute the variance reduction factor in terms of ELSS
  elss <- n_ell * var_tau_ell / var_est
  list(var = var_est, elss = elss)
}
