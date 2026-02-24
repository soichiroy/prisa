#' Prediction-error Robust Inference for linear regression
#' 
#' @param main_formula A formula for the main model. The same argument as in `lm()`.
#' @param proxy_formula A formula for the proxy model. The same argument as in `lm()`.
#' @inheritParams prisa
#' @seealso [stats::lm()]
#' @export
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
  .CallPrisaGlm(
    main_formula = main_formula,
    proxy_formula = proxy_formula,
    data = data,
    labeled_set_var_name = labeled_set_var_name,
    family = gaussian(),
    options = options,
    args_main_model = args_main_model,
    args_proxy_model = args_proxy_model
  )
}

#' Prediction-error Robust Inference for generalized linear models
#' 
#' @param main_formula A formula for the main model. The same argument as in `glm()`.
#' @param proxy_formula A formula for the proxy model. The same argument as in `glm()`.
#' @param family A description of the error distribution and link function to be
#'   used in the model. See the family argument in the `glm()` function for the 
#'   details.
#' @inheritParams prisa
#' @seealso [stats::glm()]
#' @export
prisa_glm <- function(
  main_formula,
  proxy_formula,
  data,
  labeled_set_var_name,
  family,
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
  .CallPrisaGlm(
    main_formula = main_formula,
    proxy_formula = proxy_formula,
    data = data,
    labeled_set_var_name = labeled_set_var_name,
    family = family,
    options = options,
    args_main_model = args_main_model,
    args_proxy_model = args_proxy_model
  )
}

# Internal implementation for both prisa_lm and prisa_glm
.CallPrisaGlm <- function(
  main_formula,
  proxy_formula,
  data,
  labeled_set_var_name,
  family,
  options,
  args_main_model = list(),
  args_proxy_model = list()
) {
  if (!is.list(args_main_model)) {
    stop("args_main_model must be a list.")
  }
  if (!is.list(args_proxy_model)) {
    stop("args_proxy_model must be a list.")
  }

  .main_glm_model <- function(df, formula, family, weights = NULL) {
    fit <- glm(formula, data = df, family = family, weights = weights)
    return(fit$coefficients)
  }

  .proxy_glm_model <- function(df, formula, family, weights = NULL) {
    fit <- glm(formula, data = df, family = family, weights = weights)
    return(fit$coefficients)
  }

  # Explicitly pass formula and family objects
  args_main_model$formula <- as.formula(main_formula)
  args_proxy_model$formula <- as.formula(proxy_formula)
  args_main_model$family <- family
  args_proxy_model$family <- family

  # Remove weights if present (handled elsewhere)
  if ("weights" %in% names(args_main_model)) {
    args_main_model$weights <- NULL
  }
  if ("weights" %in% names(args_proxy_model)) {
    args_proxy_model$weights <- NULL
  }

  .PredictionErrorRobustInference(
    main_model = .main_glm_model,
    proxy_model = .proxy_glm_model,
    data = data,
    labeled_set_var_name = labeled_set_var_name,
    options = options,
    args_main_model = args_main_model,
    args_proxy_model = args_proxy_model
  )
}
