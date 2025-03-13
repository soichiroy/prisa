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
#'   rows that contains the labeled data must be indicated by a variable in the
#'   data frame. The variable name must be provided in the labeled_set_var_name
#'   argument.
#' @param labeled_set_var_name The name of the variable that indicates whether 
#'  rows are labeled or not. The variable must be binary.
#' @param n_boot The number of bootstrap samples. Default is 100.
#' @param use_full A logical value that indicates whether the full data should
#'  be used in the proxy model. Default is TRUE. If FALSE, the unlabeled data
#'  will be used in the proxy model.
#' @example examples/example-MLcovar.R
#' @export
MLcovar <- function(
  main_model,
  proxy_model,
  data,
  labeled_set_var_name,
  n_boot = 100,
  use_full = TRUE
) {

  # Set options 
  # option_params <- .SetOptions(n_boot)

  # Split data into labeled and unlabeled sets
  data_list <- .SplitData(data, labeled_set_var_name)

  # Get combined estimate
  point_estimate <- .GetPointEstimates(
    main_model, proxy_model, data_list, use_full
  )

  # Naive bootstrap implementation
  boot_estimates <- .RunBootstrap(
    main_model, proxy_model, data_list, n_boot, length(point_estimate)
  )

  # Estimate optimal coefficients and variance
  coef_estimates <- .EstimateOptimalCoefficients(boot_estimates)

  # Combine estimates
  main_estimate <- .CombineEstimates(point_estimate, coef_estimates$coef)

  # Prepare output
  output <- .FormatOutput(
    data_list, main_estimate, point_estimate, coef_estimates
  )

  class(output) <- c(class(output), "MLcovar")
  return(output)
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
  return(estimates)
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

  return(list(
    dat_labeled = dat_labeled,
    dat_unlabeled = unlabeled_set,
    dat_full = data,
    n_ell = nrow(dat_labeled),
    n_full = nrow(data)
  ))

}

#' Sampling function for bootstrap
#' 
#' @param data_list A list that contains the labeled and unlabeled data.
#' @noRd 
.ConditionalSampling <- function(data_list) {
  # Conditional bootstrap to guarantee that the ratio of the labeled and
  # unlabeled sets is preserved in the resampled data. 
  #
  # TODO: Support clustered bootstrap
  
  # Sample labeled set with replacement
  dat_labeled_resampled <-
    slice_sample(data_list$dat_labeled, prop = 1, replace = TRUE)

  # Sample unlabeled set with replacement 
  dat_unlabeled_resampled <-
    slice_sample(data_list$dat_unlabeled, prop = 1, replace = TRUE)

  # Combine resampled data 
  dat_full_resampled <- 
    bind_rows(dat_labeled_resampled, dat_unlabeled_resampled)
  
  return(list(
    dat_labeled = dat_labeled_resampled,
    dat_unlabeled = dat_unlabeled_resampled,
    dat_full = dat_full_resampled
  ))
}

.RunBootstrap <- function(
    main_model, proxy_model, data_list, n_boot, n_estimates) {
  # TODO: implement parallel processing
  boot_estimate <- matrix(NA, nrow = n_boot, ncol = n_estimates)
  for (i in seq_len(n_boot)) {
    data_list_resampled <- .ConditionalSampling(data_list)
    point_estimate_resampled <-
      .GetPointEstimates(main_model, proxy_model, data_list_resampled)
    boot_estimate[i, ] <- point_estimate_resampled
  }
  return(boot_estimate)
}

#' Estimate optimal coefficients
#' 
#' @param estimates A matrix of estimates from the bootstrap samples. Each row
#'  corresponds to a bootstrap sample and each column corresponds to an estimate
#'  of the target parameter or the difference in the proxy model estimates.
#' @noRd
.EstimateOptimalCoefficients <- function(estimates) {
  # Estimate variance covariance matrix of the estimates
  cov_estimates <- cov(estimates)

  n_elements <- ncol(estimates)
  if (n_elements == 2) {
    # When the proxy estimator is also a scalar, the optimal coefficient is
    # given by Cov(tau_ell, delta_diff) / Var(delta_diff).
    coef_estimates <- cov_estimates[1, 2] / cov_estimates[2, 2]

    # Estimate variance
    var_est <- .EstimateVariance(cov_estimates) 
    return(list(coef = coef_estimates, var = var_est, vcov = cov_estimates))
  }

  # Variance of the proxy estimators
  # - This operation assumes the order of estimates in the object.
  var_diff <- cov_estimates[2:n_elements, 2:n_elements]

  # Covariance between the main and proxy model estimates
  cov_main_proxy <- as.vector(cov_estimates[1, -1])

  # Optimal coefficients
  coef_estimates <- as.vector(solve(var_diff, cov_main_proxy))

  # This check can be deleted
  if (length(coef_estimates) != (n_elements - 1)) {
    stop(
      "The number of coefficients does not match the number of proxy estimates."
    )
  }

  # Estimate variance
  var_est <- .EstimateVariance(cov_estimates) 

  return(list(coef = coef_estimates, var = var_est, vcov = cov_estimates))
}

.EstimateVariance <- function(vcov_estimate) {

  # Variance of the tau_ell
  var_tau_ell <- vcov_estimate[1, 1]

  if (ncol(vcov_estimate) == 2) {
    # Variance of the combined estimator
    var_main <- var_tau_ell - (vcov_estimate[1, 2])^2 / vcov_estimate[2, 2]
    return(var_main)
  }

  var_reduction_part <- as.vector(
    vcov_estimate[1, -1] %*% solve(vcov_estimate[-1, -1], vcov_estimate[1, -1])
  )
  return(var_tau_ell - var_reduction_part)
}

.CombineEstimates <- function(point_estimate, coef_estimates) {
  tau_ell <- point_estimate[1]
  delta_diff <- point_estimate[-1]

  # Proposed estimator: unbiased + coef * (cv_estimators)
  combined_estimate <- tau_ell + as.vector(coef_estimates %*% delta_diff) 
  return(combined_estimate)
}

.FormatOutput <- function(
    data_list, main_estimate, point_estimate, coef_estimates) {

  # Standard error of the main estimate
  std_error <- sqrt(coef_estimates$var)

  # Main information
  main_df <- data.frame(
    estimate = main_estimate,
    std_err = std_error,
    ci_lower_95 = main_estimate - qnorm(1 - 0.05 / 2) * std_error,
    ci_upper_95 = main_estimate + qnorm(1 - 0.05 / 2) * std_error,
    n_elss = .ComputeELSS(data_list$n_ell, coef_estimates)
  )

  # Other information necessary for the follow-up analysis
  output_quantities <- list(
    main_estimate = main_estimate,
    point_estimate = point_estimate,
    coef_estimates = coef_estimates
  )
  
  return(list(estimates = main_df, additional_info = output_quantities))
}