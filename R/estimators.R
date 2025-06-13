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
  tau_ell <- .FitModel(main_model, data_list$dat_labeled, args_main_model)

  # Biased estimators
  delta_ell <- .FitModel(proxy_model, data_list$dat_labeled, args_proxy_model)

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
  delta_full <- .FitModel(proxy_model, data_list$dat_full, args_proxy_model)

  if (length(delta_full) != n_proxy_estimates) {
    stop(
      "The length of the output from proxy_model must be same when evaluated on
       the labeled and full data.  Different lengths could happen for example
       when estimates include fixed effects. Consider limiting the output from
       proxy_model to main parameters of interest."
    )
  }

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

#' Fit a model using the provided function and arguments
#' 
#' @noRd
.FitModel <- function(model_function, data, args) {
  # Fit the model using the provided function and arguments
  model_fit <- do.call(model_function, c(list(data), args))
  
  # Check if the model fit is numeric and has the expected length
  if (
    !is.numeric(model_fit) || (!is.vector(model_fit) && length(model_fit) != 1)
  ) {
    stop("main_model and proxy_model functions must return a numeric scalar or vector.")
  }

  return(model_fit)
}


#' Add a name to the vector if it is not already named
#' @param x A numeric vector.
#' @noRd 
.AddName <- function(x) {
  if (!is.null(names(x)) && all(names(x) != "")) {
    return(x)
  }
  
  names(x) <- paste0("Estimate_", seq_along(x))
  return(x)
}

#' Remove names from a vector
#' @param x A numeric vector.
#' @noRd 
.RemoveName <- function(x) {
  if (is.null(names(x)) || all(names(x) == "")) {
    return(x)
  }
  
  names(x) <- NULL
  return(x)
}

#' Estimate optimal coefficients
#' 
#' @param cov_estimates An output from .ProcessCovarianceEstimates.
#' @noRd
.EstimateOptimalCoefficients <- function(cov_estimates) {
  #
  # When the vcov_full is available, estimate the coefficient using the 
  # full formula: A* = Cov(tau_main_ell, tau_proxy_diff) / V(tau_proxy_diff)
  #
  # When the vcov_full is not available, use the labeled set estimator to
  # estimate the coefficients.
  #  A* = COV(tau_proxy_ell, tau_main_ell) / V(tau_proxy_ell)
  coef_estimates <- solve(cov_estimates$vcov_delta, cov_estimates$cov_delta_tau)
  coef_estimates
}

.CombineEstimates <- function(point_estimate, coef_estimates) {
  tau_ell <- point_estimate$tau_ell
  delta_diff <- point_estimate$delta_diff


  # Proposed estimator: unbiased + coef * (cv_estimators)
  est_adjust_term <- .RemoveName(as.vector(t(coef_estimates) %*% delta_diff))
  combined_estimate <- tau_ell - est_adjust_term
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

  # Variance covariance of the main and proxy estimators
  vcov_labeled <- cov_estimates$vcov_labeled

  # Variance of the labeled-only estimator
  var_tau_ell <- diag(cov_estimates$vcov_tau) 

  # Scaling adjustment factor
  #
  # Variance of the proposed estimator:
  #   Additional (1 - prop) scaling because cov_main_proxy is based on the
  #   labeled data.
  #
  # Use the full formula to estimate the variance when use_full = TRUE.
  scale_const <- ifelse(options$use_full, 1, 1 - prop)
  vcov_reduction <- t(cov_estimates$cov_delta_tau) %*% coef_estimates
  var_est <- var_tau_ell - scale_const * diag(vcov_reduction)

  # Check the variance estimate
  var_theoretical_limit <- prop * var_tau_ell
  if (any(var_est < var_theoretical_limit) && !options$debug_mode) {
    stop(
      "Failed variance estimation. It is possible that the number of observations
       in the labeled data is too small. Please rerun with use_full = TRUE option,
       or increase the number of bootstrap iterations by increasing n_boot.
       See SetOptions for more details."
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
  list(var = var_est, elss = elss, var_labeled_only = var_tau_ell)
}
