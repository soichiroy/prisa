# Generate data
set.seed(1234)
n_ell <- 50
n_u <- 5000
n_all <- n_ell + n_u

# true model
X <- rnorm(n_all)
Y <- 0.5 + X + rnorm(n_all, mean = 0, sd = 0.3)

# Add noise to the outcome to create proxy variable
# Y_proxy are Y with non-classical measurement errors
Y_proxy_1 <- Y + rnorm(n_all, mean = X / X^2, sd = 0.1)
Y_proxy_2 <- as.numeric(Y > 0)

cor(Y, Y_proxy_1)
cor(Y, Y_proxy_2)

# Create a data frame
df_test <- data.frame(
  Y = Y,
  Y_proxy_1 = Y_proxy_1,
  Y_proxy_2 = Y_proxy_2,
  X = X,
  is_labeled = c(rep(1, n_ell), rep(0, n_u))
) |> 
  dplyr::mutate(Y = dplyr::if_else(is_labeled == 1, Y, NA_real_))

# Prepare functions for the analysis
fn_main <- function(df) {
  fit <- lm(Y ~ X, data = df)
  return(fit$coef[2])
}

fn_proxy_1 <- function(df) {
  fit <- lm(Y_proxy_1 ~ X, data = df)
  return(fit$coef)
}

fn_proxy_2 <- function(df) {
  fit <- lm(Y_proxy_2 ~ X, data = df)
  return(fit$coef)
}

fn_proxy_1_and_2 <- function(df) {
  fit1 <- lm(Y_proxy_1 ~ X, data = df)
  fit2 <- lm(Y_proxy_2 ~ X, data = df)
  return(c(fit1$coef, fit2$coef))
}

# Run functions
fit_1 <- MLcovar(
  main_model = fn_main,
  proxy_model = fn_proxy_1,
  data = df_test,
  labeled_set_var_name = "is_labeled",
  n_boot = 100
)

fit_2 <- MLcovar(
  main_model = fn_main,
  proxy_model = fn_proxy_2,
  data = df_test,
  labeled_set_var_name = "is_labeled",
  n_boot = 100
)

fit_12 <- MLcovar(
  main_model = fn_main,
  proxy_model = fn_proxy_1_and_2,
  data = df_test,
  labeled_set_var_name = "is_labeled",
  n_boot = 100
)

# Result with a noisy proxy
summary(fit_1)
sqrt(diag(fit_1$additional_info$coef_estimates$vcov))[1]

# Result with a good proxy
summary(fit_2)
sqrt(diag(fit_2$additional_info$coef_estimates$vcov))[1]

# Results with two proxies
summary(fit_12)
sqrt(diag(fit_12$additional_info$coef_estimates$vcov))[1]

# Compare with lm outputs
summary(lm(Y ~ X, data = df_test))$coef[2, 1:2]

# Biased estimates
summary(lm(Y_proxy_1 ~ X, data = df_test))$coef[2, 1:2]
summary(lm(Y_proxy_2 ~ X, data = df_test))$coef[2, 1:2]
