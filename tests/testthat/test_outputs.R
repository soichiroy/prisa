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

