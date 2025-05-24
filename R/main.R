#' @title Prediction-error Robust Inference for Statistical Analysis
#'
#' @param main_model A function that estimates the target parameter with the 
#'   labeled data. The function takes the data as the first argument and 
#'   returns the estimate of the target parameter. The target parameter must 
#'   be a scalar. The additional arguments to the main_model function must be 
#'   provided in the args_main_model argument.
#' @param proxy_model A function that estimates the target parameter with the 
#'   proxy variables. The function takes the data as the first argument and 
#'   returns the estimate of the target parameter. The target parameter must 
#'   be a scalar or a vector. The additional arguments to the proxy_model 
#'   function must be provided in the args_proxy_model argument.
#' @param data A data frame that contains the labeled and unlabeled data. The 
#'   rows that contains the labeled data must be indicated by a variable in 
#'   the data frame. The variable name must be provided in the 
#'   labeled_set_var_name argument.
#' @param labeled_set_var_name The name of the variable that indicates whether 
#'   rows are labeled or not. The variable must be binary.
#' @param options A list of options for the analysis. The values must be set 
#'   by the [SetOptions()].
#' @param args_main_model A list of additional arguments to be passed to the 
#'   main_model function. The list must be named.
#' @param args_proxy_model A list of additional arguments to be passed to the 
#'   proxy_model function. The list must be named.
#' @seealso [SetOptions()]
#' @example inst/examples/example-peri-minimal.R
#' @examplesIf interactive()
#'   # See full examples in inst/examples/example-peri.R
#'   # Or run: source(system.file("examples/example-peri.R", package = "peri"))
#' @export
peri <- function(
  main_model,
  proxy_model,
  data,
  labeled_set_var_name,
  options = SetOptions(),
  args_main_model = list(),
  args_proxy_model = list()
) {
  .PredictionErrorRobustInference(
    main_model = main_model,
    proxy_model = proxy_model,
    data = data,
    labeled_set_var_name = labeled_set_var_name,
    options = options,
    args_main_model = args_main_model,
    args_proxy_model = args_proxy_model
  )
}

#' @title Main implementation of the Prediction-error Robust Inference 
#' 
#' @inheritParams peri
#' @noRd 
.PredictionErrorRobustInference <- function(
  main_model,
  proxy_model,
  data,
  labeled_set_var_name,
  options = SetOptions(),
  args_main_model = list(),
  args_proxy_model = list()
) {

  # Split data into labeled and unlabeled sets
  data_list <- .SplitData(data, labeled_set_var_name)

  # Input check on the labeled data
  .CheckInput(data_list$dat_labeled, options)

  # Get combined estimate
  point_estimate <- .GetPointEstimates(
    main_model,
    proxy_model,
    data_list,
    args_main_model,
    args_proxy_model
  )

  # Naive bootstrap implementation
  cov_estimates <- .RunBootstrap(
    main_model,
    proxy_model,
    data_list,
    point_estimate,
    options,
    args_main_model,
    args_proxy_model
  )

  # Estimate optimal coefficients
  coef_estimates <- .EstimateOptimalCoefficients(cov_estimates)

  # Estimate the variance
  var_estimates <- .EstimateVariance(
    cov_estimates,
    coef_estimates,
    point_estimate$n_main_estimates,
    data_list$prop,
    data_list$n_ell,
    options
  )

  # Combine estimates
  main_estimate <- .CombineEstimates(point_estimate, coef_estimates)

  # Prepare output
  output <- .FormatOutput(
    main_estimate,
    point_estimate,
    coef_estimates,
    cov_estimates,
    var_estimates,
    data_list,
    options
  )

  class(output) <- c(class(output), "peri")
  output
}


#' Set options for the main function
#' 
#' Returns a list of configuration options to be passed to [peri()].
#' 
#' @param n_boot The number of bootstrap samples for the labeled data. Default 
#'  is 500. When use_full is TRUE, the full data is also used for the bootstrap.
#' @param use_full A logical value that indicates whether the full data should
#'  be used in the proxy model. Default is TRUE. If FALSE, the unlabeled data
#'  will be used in the proxy model.
#' @param use_parallel A logical value that indicates whether the bootstrap
#'  should be run in parallel. Default is TRUE. 
#' @param n_cores The number of cores to be used for bootstrap interactions.
#'  Default is the number of cores available minus one. This value will be
#'  ignored if use_parallel is FALSE.
#' @param seed_value The seed value for the random number generator. Default is
#'  drawn from a uniform between 1 and 1e7. When use_parallel is FALSE, the seed
#'  value specified here does not affect the results. 
#' @param cluster_var_name The name of the variable that indicates the cluster. 
#'   When provided, the cluster bootstrap will be used. Supports multi-way 
#'   clustering by providing a vector of variable names (e.g., c("x", "y")).
#' @param debug_mode A boolean value. If TRUE, the debug mode will be used. 
#' @return A named list of options.
#' @seealso [peri()]
#' @importFrom parallel detectCores
#' @export 
SetOptions <- function(
    n_boot = 500,
    use_full = TRUE,
    use_parallel = TRUE,
    n_cores = parallel::detectCores() - 1,
    seed_value = floor(runif(1, 1, 1e7)),
    cluster_var_name = NULL,
    debug_mode = FALSE) {

  # Check inputs
  if (!is.numeric(n_boot) || n_boot <= 0) {
    stop("n_boot must be a positive number.")
  }
  if (!is.logical(use_full)) {
    stop("use_full must be a logical value.")
  }
  if (!is.logical(use_parallel)) {
    stop("use_parallel must be a logical value.")
  }
  if (n_cores > parallel::detectCores() && use_parallel) {
    stop("n_cores must be less than the total number of cores available.")
  }
  if (!is.null(cluster_var_name) && !is.character(cluster_var_name)) {
    stop("cluster_var_name must be a string or a vector of strings.")
  }
  if (!is.null(cluster_var_name) && any(!nzchar(cluster_var_name))) {
    stop("cluster_var_name must not contain empty strings.")
  }
  if (isTRUE(debug_mode)) message("Running under the debug mode.")

  # Return the list of options
  list(
    n_boot = n_boot,
    use_full = use_full,
    use_parallel = use_parallel,
    n_cores = n_cores,
    seed_value = seed_value,
    cluster_var_name = cluster_var_name,
    debug_mode = debug_mode
  )
}

#' Split data into labeled and unlabeled sets
#' @noRd
.SplitData <- function(data, labeled_set_var_name) {
  # Check if labeled_set_var_name is present in the data
  if (!labeled_set_var_name %in% names(data)) {
    stop(paste("The variable", labeled_set_var_name, "must be in the data."))
  }
  # Check if labeled_set_var_name is binary (0 or 1)
  if (!all(data[[labeled_set_var_name]] %in% c(0, 1))) {
    stop("The labeled_set_var_name variable must be binary (0 or 1).")
  }

  # Split data into labeled and unlabeled sets
  dat_labeled <- dplyr::filter(data, !!sym(labeled_set_var_name) == 1)
  unlabeled_set <- dplyr::filter(data, !!sym(labeled_set_var_name) == 0)

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
    elss = var_estimates$elss
  )

  label_only_df <- data.frame(
    estimate = point_estimate$tau_ell,
    std_err = se_labeled_only,
    ci_lower_95 = point_estimate$tau_ell - qnorm(1 - 0.05 / 2) * se_labeled_only,
    ci_upper_95 = point_estimate$tau_ell + qnorm(1 - 0.05 / 2) * se_labeled_only
  )

  # Other information necessary for the follow-up analysis
  output_quantities <- list(
    point_estimate = point_estimate,
    coef_estimates = coef_estimates,
    cov_estimates = cov_estimates,
    var_estimates = var_estimates
  )
  
  list(
    estimates = list(main = main_df, labeled_only = label_only_df),
    additional_info = output_quantities,
    data_list = data_list[c("n_ell", "n_full", "prop")],
    options = options
  )
}
