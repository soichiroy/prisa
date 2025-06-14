
#' Implement the Prediction-error Robust Inference for linear regression
#' 
#' @param main_formula A formula for the main model. The same argument as in `lm()`.
#' @param proxy_formula A formula for the proxy model. The same argument as in `lm()`.
#' @inheritParams prisa
prisa_lm <- function(
  main_formula,
  proxy_formula,
  data,
  labeled_set_var_name,
  options = SetOptions(
    n_boot = 500,
    use_full = TRUE,
    use_parallel = TRUE,
    n_cores = parallel::detectCores() - 1,
    cluster_var_name = NULL
  ),
  args_main_model = list(),
  args_proxy_model = list()
) {
  if (!is.list(args_main_model)) {
    stop("args_main_model must be a list.")
  }
  if (!is.list(args_proxy_model)) {
    stop("args_proxy_model must be a list.")
  }
  
  # TODO: Support for additional arguments in the main and proxy models
  # TODO: Check the scope of main_formula and proxy_formula objects
  main_lm_model <- function(df, formula, weights = NULL) {
    fit <- lm(formula, data = df, weights = weights)
    return(fit$coefficients)
  }

  proxy_lm_model <- function(df, formula, weights = NULL) {
    fit <- lm(formula, data = df, weights = weights)
    return(fit$coefficients)
  }

  # Explicitly pass formula objects
  args_main_model$formula <- as.formula(main_formula)
  args_proxy_model$formula <- as.formula(proxy_formula)

  # TODO: When weights are provided, update ELSS to reflect the weights for the 
  # labeled only estimator.
  if ("weights" %in% names(args_main_model)) {
    args_main_model$weights <- NULL
  }
  if ("weights" %in% names(args_proxy_model)) {
    args_proxy_model$weights <- NULL
  }

  .PredictionErrorRobustInference(
    main_model = main_lm_model,
    proxy_model = proxy_lm_model,
    data = data,
    labeled_set_var_name = labeled_set_var_name,
    options = options,
    args_main_model = args_main_model,
    args_proxy_model = args_proxy_model
  )
}
