# Minimal example for peri function
set.seed(123)
df <- data.frame(Y = rnorm(10), Y_proxy = rnorm(10), is_labeled = c(rep(1, 5), rep(0, 5)))
fn_true <- function(df) mean(df$Y)
fn_proxy <- function(df) mean(df$Y_proxy)
peri(main_model = fn_true, proxy_model = fn_proxy, data = df, labeled_set_var_name = "is_labeled")
