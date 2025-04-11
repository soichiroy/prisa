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
#' @param options A list of options for the analysis. The values must be set by
#'  the [SetOptions()].
#' @seealso [SetOptions()]
#' @example examples/example-MLcovar.R
#' @export
MLcovar <- function(
  main_model,
  proxy_model,
  data,
  labeled_set_var_name,
  options = SetOptions()
) {

  seed_value <- options$seed_value

  # Split data into labeled and unlabeled sets
  data_list <- .SplitData(data, labeled_set_var_name)

  # Get combined estimate
  point_estimate <- .GetPointEstimates(main_model, proxy_model, data_list)

  # Naive bootstrap implementation
  cov_estimates <- .RunBootstrap(
    main_model,
    proxy_model,
    data_list,
    point_estimate,
    options
  )

  # Estimate optimal coefficients
  coef_estimates <- .EstimateOptimalCoefficients(
    cov_estimates, data_list$prop, data_list$n_ell
  )

  # Combine estimates
  main_estimate <- .CombineEstimates(point_estimate$estimates, coef_estimates)

  # Estimate the variance
  var_estimates <- .EstimateVariance(
    cov_estimates, coef_estimates, data_list$prop, data_list$n_ell
  )

  # Prepare output
  output <- .FormatOutput(
    main_estimate,
    point_estimate$estimates,
    coef_estimates,
    cov_estimates,
    var_estimates,
    data_list,
    options
  )

  class(output) <- c(class(output), "MLcovar")
  output
}


#' Set options for the main function
#' 
#' Returns a list of configuration options to be passed to [MLcovar()].
#' 
#' @param n_boot The number of bootstrap samples for the labeled data. Default 
#'  is 500. When use_full is TRUE, the full data is also used for the bootstrap.
#' @param use_full A logical value that indicates whether the full data should
#'  be used in the proxy model. Default is TRUE. If FALSE, the unlabeled data
#'  will be used in the proxy model.
#' @param use_parallel A logical value that indicates whether the bootstrap
#'  should be run in parallel. Default is TRUE. 
#' @param seed_value The seed value for the random number generator. Default is
#'  drawn from a uniform between 1 and 1e7. 
#' 
#' @return A named list of options.
#' @seealso [MLcovar()]
#' @export 
SetOptions <- function(
    n_boot = 500,
    use_full = TRUE,
    use_parallel = TRUE,
    seed_value = floor(runif(1, 1, 1e7))) {
  list(
    n_boot = n_boot,
    use_full = use_full,
    use_parallel = use_parallel,
    seed_value = seed_value
  )
}

.GetPointEstimates <- function(main_model, proxy_model, data_list) {
  # Unbiased estimator
  tau_ell <- main_model(data_list$dat_labeled)

  if (length(as.vector(tau_ell)) > 1) {
    stop("The main_model function must return a scalar value.")
  }

  # Biased estimators
  delta_ell <- proxy_model(data_list$dat_labeled)
  delta_full <- proxy_model(data_list$dat_full)

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

.GetPointEstimatesLabeled <- function(
    main_model, proxy_model, data_labeled_resampled
) {

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

  n_ell <- nrow(dat_labeled)
  n_full <- nrow(data)
  prop <- n_ell / n_full
  list(
    dat_labeled = dat_labeled,
    dat_unlabeled = unlabeled_set,
    dat_full = data,
    n_ell = n_ell,
    n_full = n_full,
    prop = prop
  )
}

#' Run bootstrap
#' 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doRNG registerDoRNG
#' @import dplyr
#' @noRd 
.RunBootstrap <- function(
    main_model,
    proxy_model,
    data_list,
    point_estimate,
    options) {
  
  # Extract option values
  n_boot <- options$n_boot
  use_full <- options$use_full
  use_parallel <- options$use_parallel
  seed_value <- options$seed_value

  n_estimates_labeled <- point_estimate$n_estimates_labeled
  n_estimates_full <- point_estimate$n_estimates_full

  # Register parallel backend
  if (use_parallel) {
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    doRNG::registerDoRNG(seed_value)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    })
  } else {
    # When use_parallel is FALSE, use sequential processing
    foreach::registerDoSEQ()
  }

  # Bootstrap for the labeled data to estimate the variance-covariance matrix
  # of the main estimators based on the labeled data.
  boot_estimate_labeled <- foreach(
    i = seq_len(n_boot),
    .combine = rbind,
    .packages = c("dplyr"),
    .inorder = FALSE,
    .export = c(".GetPointEstimatesLabeled")
  ) %dopar% {
      data_labeled_resampled <- 
        slice_sample(data_list$dat_labeled, prop = 1, replace = TRUE)
      .GetPointEstimatesLabeled(
        main_model, proxy_model, data_labeled_resampled
      )
  }
  vcov_labeled <- cov(boot_estimate_labeled)

  # Exit the function if boot_full is FALSE (skip the bootstrap for the entire
  # data).
  if (isFALSE(use_full)) {
    # Return vcov estimate
    return(list(
      vcov_labeled = vcov_labeled,
      vcov_full = NULL,
      vcov_main_diff = NULL
    ))
  }


  # Estimate the variance covariance matrix of the biased estimator based on 
  # the unlabeled data or the full data
  data_main <- data_list$dat_full
  
  # Bootstrap for the full (or unlabeled) data
  boot_estimate_full <- foreach(
    i = seq_len(n_boot),
    .combine = rbind,
    .inorder = FALSE,
    .packages = c("dplyr")
  ) %dopar% {
      data_main_resampled <- slice_sample(data_main, prop = 1, replace = TRUE)
      proxy_model(data_main_resampled)
  }

  # VCOV(tau_proxy_ell - tau_proxy_full)
  boot_estimate_diff <- boot_estimate_labeled[, -1] - boot_estimate_full
  vcov_full <- cov(boot_estimate_diff)

  # Cov(tau_main_ell, tau_proxy_ell - tau_proxy_full)
  vcov_main_diff <- 
    as.vector(cov(boot_estimate_labeled[, 1], boot_estimate_diff))

  # Return bootstrapped variance-covariance estimates
  list(
    vcov_labeled = vcov_labeled,
    vcov_full = vcov_full,
    vcov_main_diff = vcov_main_diff
  )
}

#' Estimate optimal coefficients
#' 
#' @param cov_estimates A list of bootstrap variance-covariance estimates.
#' @noRd
.EstimateOptimalCoefficients <- function(cov_estimates, prop, n_ell) {
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
.EstimateVariance <- function(cov_estimates, coef_estimates, prop, n_ell) {
  # Variance covariance of the main and proxy estimators
  vcov_labeled <- cov_estimates$vcov_labeled
  cov_main_proxy <- as.vector(vcov_labeled[1, -1])

  # Variance of the original estimator
  var_tau_ell <- vcov_labeled[1, 1]

  if (is.null(cov_estimates$vcov_full)) {
    # Variance of the proposed estimator:
    #   Additional (1 - prop) scaling because cov_main_proxy is based on the
    #   labeled data.
    var_est <- var_tau_ell - 
      (1 - prop) * as.vector(cov_main_proxy %*% coef_estimates) 
  } else {
    var_est <- var_tau_ell - 
      as.vector(cov_estimates$vcov_main_diff %*% coef_estimates)
  }

  var_theoretical_limit <- prop * var_tau_ell
  if (var_est < var_theoretical_limit) {
    warning("Estimated variance is too small. Using the theoretical limit.")
  	var_est <- var_theoretical_limit
  }

  # Compute the variance reduction factor in terms of ELSS
  elss <- n_ell * var_tau_ell / var_est
  list(var = var_est, elss = elss)
}

#' Format the output of the proposed estimator
#' 
#' @noRd
.FormatOutput <- function(
    main_estimate,
    point_estimate,
    coef_estimates,
    cov_estimates,
    var_estimates,
    data_list,
    options) {

  # Standard error of the proposed estimator
  std_error <- sqrt(var_estimates$var)

  # Standard error of the labeled-only estimator
  se_labeled_only <- sqrt(cov_estimates$vcov_labeled[1, 1])

  # Main information
  main_df <- data.frame(
    estimate = main_estimate,
    std_err = std_error,
    ci_lower_95 = main_estimate - qnorm(1 - 0.05 / 2) * std_error,
    ci_upper_95 = main_estimate + qnorm(1 - 0.05 / 2) * std_error,
    elss = var_estimates$elss,
    labeled_only = point_estimate[1],
    se_labeled_only = se_labeled_only
  )

  # Other information necessary for the follow-up analysis
  output_quantities <- list(
    point_estimate = point_estimate,
    coef_estimates = coef_estimates,
    cov_estimates = cov_estimates,
    var_estimates = var_estimates
  )
  
  list(
    estimates = main_df,
    additional_info = output_quantities,
    data_list = data_list[c("n_ell", "n_full", "prop")],
    options = options
  )
}