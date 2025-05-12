df_list <- list(
  dat_labeled = data.frame(x = 1:5, y = (1:5) / 2),
  dat_full = data.frame(y = (1:40) / 2)
)
main_model <- function(df) {
  return(mean(df$x))
}
proxy_model <- function(df) {
  return(median(df$y))
}

test_that("Return correct output from .GetPointEstimates", {
  expect_equal(.GetPointEstimates(
    main_model = main_model,
    proxy_model = proxy_model,
    data_list = df_list,
    args_main_model = list(),
    args_proxy_model = list()
  ), list(
    estimates = c(mean(1:5), (mean(1:5) - mean(1:40)) / 2),
    # labeled-only estimates + proxy estimates
    n_estimates_labeled = 2,
    # proxy estimates
    n_estimates_full = 1
  ))
})

test_that("Return correct output from .GetPointEstimatesLabeled", {
  labeled_data <- df_list$dat_labeled

  expect_equal(.GetPointEstimatesLabeled(
    main_model = main_model,
    proxy_model = proxy_model,
    data_labeled = labeled_data,
    args_main_model = list(),
    args_proxy_model = list()
  ), c(mean(1:5), mean(1:5) / 2))
})

test_that("Check valid form of functions", {
  main_model_invalid <- function(df) {
    return(c(mean(df$x), mean(df$y)))
  }
  expect_error(.GetPointEstimates(
    main_model = main_model_invalid, 
    proxy_model = proxy_model,
    data_list = df_list,
    args_main_model = list(),
    args_proxy_model = list()
  ), "The main_model function must return a scalar value.")
})

# test_that("Estimate optimal coefficients when vcov_full is available", {
#   cov_estimates <- list(
#     vcov_labeled = matrix(c(1, 0.5, 0.5, 2), nrow = 2),
#     vcov_full = matrix(c(2, 0.3, 0.3, 1), nrow = 2),
#     vcov_main_diff = c(0.4, 0.6)
#   )
#   prop <- 0.5
#   n_ell <- 100

#   coef_estimates <- .EstimateOptimalCoefficients(cov_estimates, prop, n_ell)

#   # Expected coefficients
#   expected_coef <- solve(cov_estimates$vcov_full, cov_estimates$vcov_main_diff)
#   expect_equal(coef_estimates, as.vector(expected_coef))
# })

# test_that("Estimate optimal coefficients when vcov_full is not available", {
#   cov_estimates <- list(
#     vcov_labeled = matrix(c(1, 0.5, 0.5, 2), nrow = 2),
#     vcov_full = NULL
#   )
#   prop <- 0.5
#   n_ell <- 100

#   coef_estimates <- .EstimateOptimalCoefficients(cov_estimates, prop, n_ell)

#   # Expected coefficients
#   cov_main_proxy <- as.vector(cov_estimates$vcov_labeled[1, -1])
#   expected_coef <- solve(cov_estimates$vcov_labeled[-1, -1], cov_main_proxy)
#   expect_equal(coef_estimates, as.vector(expected_coef))
# })

# test_that("Handle invalid inputs in .EstimateOptimalCoefficients", {
#   cov_estimates <- list(
#     vcov_labeled = matrix(c(1, 0.5, 0.5, 2), nrow = 2),
#     vcov_full = NULL
#   )
#   prop <- 0.5
#   n_ell <- 100

#   # Invalid vcov_labeled (non-invertible matrix)
#   cov_estimates$vcov_labeled <- matrix(c(1, 1, 1, 1), nrow = 2)
#   expect_error(
#     .EstimateOptimalCoefficients(cov_estimates, prop, n_ell),
#     "Lapack routine dgesv: system is exactly singular"
#   )
# })