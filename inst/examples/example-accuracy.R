set.seed(1234)
n_ell <- 200
n_u <- 5000
n_all <- n_ell + n_u

# true model
X <- rnorm(n_all)
Y <- 0.5 + X + rnorm(n_all, mean = 0, sd = 0.3)

# Add noise to the outcome to create proxy variable
# Y_proxy are Y with non-classical measurement errors
Y_proxy_1 <- Y + rnorm(n_all, mean = X, sd = 0.1)

# Create a data frame
df_test <- data.frame(
  Y = Y,
  Y_proxy_1 = Y_proxy_1,
  X = X,
  is_labeled = c(rep(1, n_ell), rep(0, n_u))
) |> 
  dplyr::mutate(Y = dplyr::if_else(is_labeled == 1, Y, NA_real_))

# Prepare functions for the analysis
fn_main <- function(df) {
  fit <- lm(Y ~ X, data = df)
  fit$coef[2]
}

fn_proxy_1 <- function(df) {
  fit <- lm(Y_proxy_1 ~ X, data = df)
  fit$coef
}

# Run functions
fit_1 <- MLcovar(
  main_model = fn_main,
  proxy_model = fn_proxy_1,
  data = df_test,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 500, use_full = FALSE)
)
summary(fit_1)

accuracy_out <- accuracy(fit_1)
plot(accuracy_out)