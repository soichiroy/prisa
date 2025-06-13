
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

  # TODO: Support for additional arguments in the main and proxy models
  # TODO: Check the scope of main_formula and proxy_formula objects
  main_lm_model <- function(df, formula) {
    fit <- lm(formula, data = df)
    return(fit$coefficients)
  }

  proxy_lm_model <- function(df, formula) {
    fit <- lm(formula, data = df)
    return(fit$coefficients)
  }

  # Explicitly pass formula objects
  args_main_model <-
    c(args_main_model, list(formula = as.formula(main_formula)))
  args_proxy_model <-
    c(args_proxy_model, list(formula = as.formula(proxy_formula)))

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
