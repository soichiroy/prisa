# Prepare test data
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
  expect_equal(
    .GetPointEstimates(
      main_model = main_model,
      proxy_model = proxy_model,
      data_list = df_list,
      args_main_model = list(),
      args_proxy_model = list()
    ),
    list(
      estimates = c(mean(1:5), (mean(1:5) - mean(1:40)) / 2),
      # labeled-only estimates + proxy estimates
      n_estimates_labeled = 2,
      # proxy estimates
      n_estimates_full = 1
    )
  )
})

test_that("Check valid form of functions", {
  main_model_invalid <- function(df) {
    return(c(mean(df$x), mean(df$y)))
  }
  expect_error(
    .GetPointEstimates(
      main_model = main_model_invalid,
      proxy_model = proxy_model,
      data_list = df_list,
      args_main_model = list(),
      args_proxy_model = list()
    ),
    "The main_model function must return a scalar value."
  )
})
