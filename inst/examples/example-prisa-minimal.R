# Minimal example for prisa function
set.seed(1234)

# True labels (unobserved)
Y_true <- rnorm(5100)
df <- data.frame(
  Y_label = c(Y_true[1:100], rep(NA, 5000)), 
  Y_proxy = round(Y_true / 4) * 4 + rnorm(n = 5100, mean = 0, sd = 0.5),
  is_labeled = c(rep(1, 100), rep(0, 5000))
)

# Main model is to take the sample mean of the labeled outcome
fn_true <- function(df) {
  mean(df$Y_label, na.rm = TRUE)
}

# Proxy model is to take the sample mean of the predicted outcome
fn_proxy <- function(df) {
  mean(df$Y_proxy, na.rm = TRUE)
}

fit <- prisa(
  main_model = fn_true,
  proxy_model = fn_proxy,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(
    n_boot = 100, 
    use_full = TRUE,
    use_parallel = FALSE
  )
)

summary(fit)
