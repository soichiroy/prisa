# Prepare test data
df_list <- list(
  dat_labeled = data.frame(x = 1:5, y = (1:5) / 2),
  dat_full = data.frame(y = (1:40) / 2)
)
main_model <- function(df) {
  return(mean(df$x))
}
proxy_model <- function(df) {
  return(mean(df$y))
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
      tau_ell = mean(1:5),
      delta_diff = (mean(1:5) - mean(1:40)) / 2,
      n_main_estimates = 1,
      n_proxy_estimates = 1
    )
  )
})
