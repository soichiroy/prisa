#' @title Estimators for the main and proxy models
#'
#' @noRd 
#' @inheritParams MLcovar
.GetPointEstimates <- function(
    main_model,
    proxy_model,
    data_list,
    args_main_model,
    args_proxy_model) {

  # Unbiased estimator
  tau_ell <-
    do.call(main_model, c(list(data_list$dat_labeled), args_main_model))

  if (length(as.vector(tau_ell)) > 1) {
    stop("The main_model function must return a scalar value.")
  }

  # Biased estimators
  delta_ell <- do.call(
    proxy_model, c(list(data_list$dat_labeled), args_proxy_model)
  )
  delta_full <- do.call(
    proxy_model, c(list(data_list$dat_full), args_proxy_model)
  )

  # Create control variate estimators
  delta_diff <- as.vector(delta_ell) - as.vector(delta_full)
  
  # Return estimates: Assume the first element is the main model estimate
  estimates <- c(tau_ell, delta_diff)
  n_estimates_labeled <- length(estimates)
  n_estimates_full <- length(delta_diff)
  list(
    estimates = estimates, 
    n_estimates_labeled = n_estimates_labeled, 
    n_estimates_full = n_estimates_full
  )
}


# TODO: Update .GetPointEstimates to replace this function.
.GetPointEstimatesLabeled <- function(
    main_model,
    proxy_model,
    data_labeled_resampled,
    args_main_model,
    args_proxy_model
) {

  # Unbiased estimator
  tau_ell <- do.call(
    main_model, c(list(data_labeled_resampled), args_main_model)
  )

  # Biased estimators
  delta_ell <- do.call(
    proxy_model, c(list(data_labeled_resampled), args_proxy_model)
  )

  # Return estimates
  estimates <- c(tau_ell, delta_ell) 
  estimates
}


#' Estimate optimal coefficients
#' 
#' @param cov_estimates A list of bootstrap variance-covariance estimates.
#' @noRd
.EstimateOptimalCoefficients <- function(cov_estimates) {
  vcov_labeled <- cov_estimates$vcov_labeled
  # Covariance between the main and proxy model estimates
  #  * [1,1] element is the variance of the main unbiased estimator 
  cov_main_proxy <- as.vector(vcov_labeled[1, -1])
  
  # Variance covariance matrix of the proxy estimator
  vcov_full <- cov_estimates$vcov_full

  if (!is.null(vcov_full)) {
    # When the vcov_full is available, estimate the coefficient using the 
    # following formula:
    #   A* = Cov(tau_main_ell, tau_proxy_diff) / V(tau_proxy_diff)
    coef_estimates <- as.vector(solve(vcov_full, cov_estimates$vcov_main_diff))
    return(coef_estimates)
  }

  # When the vcov_full is not available, use the labeled set estimator to
  # estimate the coefficients.
  #  A* = COV(tau_proxy_ell, tau_main_ell) / V(tau_proxy_ell)
  coef_estimates <- as.vector(solve(vcov_labeled[-1, -1], cov_main_proxy))
  coef_estimates
}

.CombineEstimates <- function(point_estimate, coef_estimates) {
  tau_ell <- point_estimate[1]
  delta_diff <- point_estimate[-1]

  # Proposed estimator: unbiased + coef * (cv_estimators)
  combined_estimate <- tau_ell - as.vector(coef_estimates %*% delta_diff) 
  combined_estimate
}

#' Estimate the variance of the proposed estimator
#' 
#' @noRd
.EstimateVariance <- function(
  cov_estimates, coef_estimates, prop, n_ell, options
) {
  # Variance covariance of the main and proxy estimators
  vcov_labeled <- cov_estimates$vcov_labeled
  cov_main_proxy <- as.vector(vcov_labeled[1, -1])

  # Variance of the original estimator
  var_tau_ell <- vcov_labeled[1, 1]

  # vcov_main_diff is NULL when use_full is FALSE.
  if (is.null(cov_estimates$vcov_main_diff)) {
    # Variance of the proposed estimator:
    #   Additional (1 - prop) scaling because cov_main_proxy is based on the
    #   labeled data.
    var_est <- var_tau_ell - 
      (1 - prop) * as.vector(cov_main_proxy %*% coef_estimates) 
  } else {
    # Use the full formula to estimate the variance
    var_est <- var_tau_ell - 
      as.vector(cov_estimates$vcov_main_diff %*% coef_estimates)
  }

  # Check the variance estimate
  var_theoretical_limit <- prop * var_tau_ell
  if (var_est < var_theoretical_limit && !options$debug_mode) {
    stop(
      "Failed variance estimation. Please use use_full = TRUE option instead."
    )
  }
   
  if (var_est < var_theoretical_limit && options$debug_mode) {
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
