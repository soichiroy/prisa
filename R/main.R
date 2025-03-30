# To be renamed after the package name is finalized


#' @title Proposed estimator
#' 
#' @param main_model A function that estimates the target parameter with the 
#'  labeled data. The function takes the data as the only argument and returns
#'  the estimate of the target parameter. The target parameter must be a scalar.
#' @param proxy_model A function that estimates the target parameter with the
#'  proxy variables. The function takes the data as the only argument and
#'  returns the estimate of the target parameter. The target parameter must be a
#'  scalar or a vector.
#' @param data A data frame that contains the labeled and unlabeled data. The 
#'  rows that contains the labeled data must be indicated by a variable in the
#'  data frame. The variable name must be provided in the labeled_set_var_name
#'  argument.
#' @param labeled_set_var_name The name of the variable that indicates whether 
#'  rows are labeled or not. The variable must be binary.
#' @param n_boot The number of bootstrap samples for the labeled data. Default 
#'  is 500.
#' @param n_boot2 The number of bootstrap samples for the full data. Default is 
#'  100.
#' @param use_full A logical value that indicates whether the full data should
#'  be used in the proxy model. Default is TRUE. If FALSE, the unlabeled data
#'  will be used in the proxy model.
#' @param boot_full A logical value that indicates whether the full data should
#'  be used in estimating the variance of the parameters with the proxy 
#'  variables for the full data. Default is TRUE. If FALSE, the unlabeled data 
#'  will be used. This option also affects point estimates.
#' @example examples/example-MLcovar.R
#' @export
MLcovar <- function(
  main_model,
  proxy_model,
  data,
  labeled_set_var_name,
  n_boot = 500,
  n_boot2 = 100,
  use_full = TRUE,
  boot_full = TRUE
) {

  # Set options 
  # option_params <- .SetOptions(n_boot)

  # Split data into labeled and unlabeled sets
  data_list <- .SplitData(data, labeled_set_var_name)

  # Proportion of the labeled data
  prop <- data_list$n_ell / data_list$n_full

  # Get combined estimate
  point_estimate <- .GetPointEstimates(
    main_model, proxy_model, data_list, use_full
  )

  # Naive bootstrap implementation
  cov_estimates <- .RunBootstrap(
    main_model,
    proxy_model,
    data_list,
    n_boot,
    n_boot2,
    point_estimate$n_estimates_labeled,
    point_estimate$n_estimates_full,
    use_full,
    boot_full
  )

  # Estimate optimal coefficients and variance
  coef_estimates <- .EstimateOptimalCoefficients(
    cov_estimates, prop, data_list$n_ell, use_full, boot_full
  )

  # Combine estimates
  main_estimate <- .CombineEstimates(
    point_estimate$estimates, coef_estimates$coef
  )

  # Prepare output
  output <- .FormatOutput(
    main_estimate,
    point_estimate$estimates,
    coef_estimates,
    cov_estimates,
    data_list
  )

  class(output) <- c(class(output), "MLcovar")
  output
}

.GetPointEstimates <- function(main_model, proxy_model, data_list, use_full) {

  # Unbiased estimator
  tau_ell <- main_model(data_list$dat_labeled)

  if (length(as.vector(tau_ell)) > 1) {
    stop("The main_model function must return a scalar value.")
  }

  # Biased estimators
  delta_ell <- proxy_model(data_list$dat_labeled)

  if (isTRUE(use_full)) {
    delta_full <- proxy_model(data_list$dat_full)
  } else {
    delta_full <- proxy_model(data_list$dat_unlabeled)
  }

  delta_diff <- as.vector(delta_ell) - as.vector(delta_full)
  
  # Return estimates
  estimates <- c(tau_ell, delta_diff)
  n_estimates_labeled <- length(estimates)
  n_estimates_full <- length(delta_diff)
  list(
  	estimates = estimates, 
  	n_estimates_labeled = n_estimates_labeled, 
  	n_estimates_full = n_estimates_full
  )
}

.GetPointEstimatesLabeled <- function(
    main_model, proxy_model, data_labeled_resampled) {

  # Unbiased estimator
  tau_ell <- main_model(data_labeled_resampled)

  # Biased estimators
  delta_ell <- proxy_model(data_labeled_resampled)

  # Return estimates
  estimates <- c(tau_ell, delta_ell) 
  estimates
}

.SplitData <- function(data, labeled_set_var_name) {

  # Check if labeled_set_var_name is binary (0 or 1)
  if (!all(data[[labeled_set_var_name]] %in% c(0, 1))) {
    stop("The labeled_set_var_name variable must be binary (0 or 1).")
  }

  # Split data into labeled and unlabeled sets
  dat_labeled <- data %>%
    filter(!!sym(labeled_set_var_name) == 1)
  unlabeled_set <- data %>%
    filter(!!sym(labeled_set_var_name) == 0)

  list(
    dat_labeled = dat_labeled,
    dat_unlabeled = unlabeled_set,
    dat_full = data,
    n_ell = nrow(dat_labeled),
    n_full = nrow(data)
  )
}

.RunBootstrap <- function(
    main_model, proxy_model, data_list, n_boot, n_boot2, n_estimates_labeled, 
    n_estimates_full, use_full, boot_full) {
  # TODO: implement parallel processing
  # Bootstrap for the labeled data
  boot_estimate_labeled <- matrix(NA, nrow = n_boot, ncol = n_estimates_labeled)
  for (i in seq_len(n_boot)) {
    data_labeled_resampled <-
      slice_sample(data_list$dat_labeled, prop = 1, replace = TRUE)
    point_estimate_boot_labeled <- .GetPointEstimatesLabeled(
      main_model, proxy_model, data_labeled_resampled
    )
    boot_estimate_labeled[i, ] <- point_estimate_boot_labeled
  }

  # Estimate variance-covariance matrix of the estimators for the labeled data
  vcov_labeled <- cov(boot_estimate_labeled)

  if (isTRUE(boot_full)) {
	  # Bootstrap for the full (or unlabeled) data
	  boot_estimate_full <- matrix(NA, nrow = n_boot2, ncol = n_estimates_full)
	  for (i in seq_len(n_boot2)) {
		  if (isTRUE(use_full)) {
		    data_full_resampled <-
		      slice_sample(data_list$dat_full, prop = 1, replace = TRUE)
		  } else {
		    data_full_resampled <-
		      slice_sample(data_list$dat_unlabeled, prop = 1, replace = TRUE)
		  }
	    boot_estimate_full[i, ] <- proxy_model(data_full_resampled)
	  }
	} else {
		vcov_full <- NULL
	}

  # Estimate variance-covariance matrix of the estimators for the full data
  vcov_full <- cov(boot_estimate_full)

  # Return bootstrapped variance-covariance estimates
  list(
  	vcov_labeled = vcov_labeled, 
  	vcov_full = vcov_full
  )
}

#' Estimate optimal coefficients
#' 
#' @param cov_estimates A list of bootstrap variance-covariance estimates.
#' @noRd
.EstimateOptimalCoefficients <- function(
    cov_estimates, prop, n_ell, use_full, boot_full) {
	vcov_labeled <- cov_estimates$vcov_labeled
	vcov_full <- cov_estimates$vcov_full
  n_elements <- ncol(cov_estimates$vcov_labeled)
	if(isFALSE(boot_full)) {
		if(isTRUE(use_full)) {
			# Estimate variance of the estimators for the full data
			vcov_full <- prop * vcov_labeled[2:n_elements, 2:n_elements]
		} else {
			# Estimate variance of the estimators for the unlabeled data
			vcov_full <- prop / (1 - prop) * vcov_labeled[2:n_elements, 2:n_elements]
		}
	}

#  if (n_elements == 2) {
#    # When the proxy estimator is also a scalar
#    coef_estimates <- prop * vcov_labeled[1, 2] / vcov_full
#
#    # Estimate variance of the combined estimator
#    var_tau_ell <- vcov_labeled[1, 1]
#    var_est <- var_tau_ell - prop * (1 - prop) * vcov_labeled[1, 2]^2 / vcov_full
#    var_theoretical_limit <- prop * var_tau_ell
#	   if (var_est < var_theoretical_limit) {
#	  	 var_est <- var_theoretical_limit
#	   }
#    elss <- n_ell * var_tau_ell / var_est
#    return(list(coef = coef_estimates, var = var_est, elss = elss))
#  }

  # Covariance between the main and proxy model estimates
  cov_main_proxy <- as.vector(vcov_labeled[1, -1])

  # Optimal coefficients
  coef_estimates <- as.vector(prop * (1 - prop) * solve(vcov_full, cov_main_proxy))

  # This check can be deleted
  if (length(coef_estimates) != (n_elements - 1)) {
    stop(
      "The number of coefficients does not match the number of proxy estimates."
    )
  }

  # Estimate variance of the combined estimator
  var_tau_ell <- vcov_labeled[1, 1]
  var_est <- var_tau_ell - 
    prop * (1 - prop) * cov_main_proxy %*% solve(vcov_full, cov_main_proxy)
  var_theoretical_limit <- prop * var_tau_ell
  if (var_est < var_theoretical_limit) {
  	var_est <- var_theoretical_limit
  }
  elss <- n_ell * var_tau_ell / var_est

  list(
  	coef = coef_estimates, 
  	var = var_est, 
  	elss = elss
  )
}

.CombineEstimates <- function(point_estimate, coef_estimates) {
  tau_ell <- point_estimate[1]
  delta_diff <- point_estimate[-1]

  # Proposed estimator: unbiased + coef * (cv_estimators)
  combined_estimate <- tau_ell - as.vector(coef_estimates %*% delta_diff) 
  combined_estimate
}

.FormatOutput <- function(
    main_estimate, point_estimate, coef_estimates, cov_estimates, data_list) {

  # Standard error of the main estimate
  std_error <- sqrt(coef_estimates$var)

  # Main information
  main_df <- data.frame(
    estimate = main_estimate,
    std_err = std_error,
    ci_lower_95 = main_estimate - qnorm(1 - 0.05 / 2) * std_error,
    ci_upper_95 = main_estimate + qnorm(1 - 0.05 / 2) * std_error,
    elss = coef_estimates$elss
  )

  # Other information necessary for the follow-up analysis
  output_quantities <- list(
    point_estimate = point_estimate,
    coef_estimates = coef_estimates,
    cov_estimates = cov_estimates
  )
  
  list(estimates = main_df, additional_info = output_quantities)
}