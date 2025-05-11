## Function
logistic <- function (x) {
  exp(x) / (1 + exp(x))
}

## Generate data
## Setup
set.seed(2025)
n_ell <- 200
n_u <- 5000
n_all <- n_ell + n_u

## Covariates
X1 <- rnorm(n = n_all, mean = 0, sd = 1)
X2 <- rnorm(n = n_all, mean = 0, sd = 1)

## Add noise to the covariates to create proxy covariates
X1_proxy <- X1 + rnorm(n = n_all, mean = 0, sd = 0.2)
X2_proxy <- X2 + rnorm(n = n_all, mean = 0, sd = 0.2)

## Instrument variable
Z <- rnorm(n = n_all, mean = 0, sd = 1)

## Add noise to the covariates to create proxy covariates
## Z_proxy does not satisfy the exclusion restriction
Z_proxy <- Z + rnorm(n = n_all, mean = (X1 + X2) / 5, sd = 0.1)

## Treatment indicator
prob_d <- logistic(-0.3 + 2 * Z + 0.4 * X1 + 0.2 * X2)
D <- rbinom(n = n_all, size = 1, prob = prob_d)

## Add noise to the treatment to create a proxy treatment (misclassifiction prob.: 3%)
D_proxy <- rbinom(n = n_all, size = 1, prob = 0.03 + D * 0.94)

## True outcome model for linear regression model
Y <- 0.5 + D + X1 + 1.5 * X2 + rnorm(n = n_all, mean = 0, sd = 1)

## Add noise to the outcome to create proxy outcomes
## Y_proxy_1 and Y_proxy_2 are Y with non-classical measurement errors
Y_proxy_1 <- Y + rnorm(n = n_all, mean = D * X1 + (1 - D) * X2, sd = 1)
Y_proxy_2 <- round(Y / 4) * 4 + rnorm(n = n_all, mean = 0, sd = 0.5)

## Check correlation
cor(Y, Y_proxy_1)
cor(Y, Y_proxy_2)

## True outcome model for a binary outcome
prob_y <- logistic(-0.5 + 2 * D + X1 + X2)
Y_binary <- rbinom(n = n_all, size = 1, prob = prob_y)

## Add noise to the binary outcome to create proxy binary outcomes (misclassifiction prob.: 5%)
Y_binary_proxy <- rbinom(n = n_all, size = 1, prob = 0.05 + Y_binary * 0.9)


# Create a data frame
df <- data.frame(
  Y = Y,
  Y_proxy_1 = Y_proxy_1,
  Y_proxy_2 = Y_proxy_2,
  Y_binary = Y_binary,
  Y_binary_proxy = Y_binary_proxy,
  D = D,
  D_proxy = D_proxy,
  Z = Z,
  Z_proxy = Z_proxy,
  X1 = X1,
  X2 = X2,
  X1_proxy = X1_proxy,
  X2_proxy = X2_proxy,
  is_labeled = c(rep(1, n_ell), rep(0, n_u))
) |> 
  dplyr::mutate(Y = dplyr::if_else(is_labeled == 1, Y, NA_real_),
                D = dplyr::if_else(is_labeled == 1, D, NA_real_),
                X1 = dplyr::if_else(is_labeled == 1, X1, NA_real_),
                X2 = dplyr::if_else(is_labeled == 1, X2, NA_real_),
                Z = dplyr::if_else(is_labeled == 1, Z, NA_real_))

## Estimation
## 0. Estimating the outcome mean
## Prepare an estimation model with true labels
fn_true_mean <- function (df) {
  mean(df$Y)
}

## Prepare an estimation model with proxy labels
fn_proxy_mean <- function (df) {
  mean(df$Y_proxy_1)
}

## Estimate with the proposed estimator
fit_mean <- MLcovar(
  main_model = fn_true_mean,
  proxy_model = fn_proxy_mean,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_mean)

## True value
mean(0.5 + D + X1 + 1.5 * X2)


## 1. Estimating the linear regression coefficient
## Prepare an estimation model with true labels
fn_true_lm <- function (df) {
  fit <- lm(Y ~ D + X1 + X2, data = df)
  fit$coefficients["D"]
}

## Prepare estimation models with proxy labels
fn_proxy_lm_1 <- function (df) {
  fit <- lm(Y_proxy_1 ~ D_proxy + X1_proxy + X2_proxy, data = df)
  fit$coefficients
}

fn_proxy_lm_2 <- function (df) {
  fit <- lm(Y_proxy_2 ~ D_proxy + X1_proxy + X2_proxy, data = df)
  fit$coefficients
}

fn_proxy_lm_12 <- function (df) {
  fit1 <- lm(Y_proxy_1 ~ D_proxy + X1_proxy + X2_proxy, data = df)
  fit2 <- lm(Y_proxy_2 ~ D_proxy + X1_proxy + X2_proxy, data = df)
  c(fit1$coefficients, fit2$coefficients)
}

fn_proxy_lm <- function(df, method = c("proxy1, proxy2", "proxy_1_2")) {
  method <- match.arg(method)
  if (method == "proxy1") {
    fit <- lm(Y_proxy_1 ~ D_proxy + X1_proxy + X2_proxy, data = df)
    return(fit$coef)
  } else if (method == "proxy2") {
    fit <- lm(Y_proxy_2 ~ D_proxy + X1_proxy + X2_proxy, data = df)
    return(fit$coef)
  }

  fit1 <- lm(Y_proxy_1 ~ D_proxy + X1_proxy + X2_proxy, data = df)
  fit2 <- lm(Y_proxy_2 ~ D_proxy + X1_proxy + X2_proxy, data = df)
  c(fit1$coefficients, fit2$coefficients)
}

## Estimate with the proposed estimator
## With a model for Y_proxy_1
fit_lm_1 <- MLcovar(
  main_model = fn_true_lm,
  proxy_model = fn_proxy_lm,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE),
  args_proxy_model = list(method = "proxy1")
)

## With a model for Y_proxy_2
fit_lm_2 <- MLcovar(
  main_model = fn_true_lm,
  proxy_model = fn_proxy_lm,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE),
  args_proxy_model = list(method = "proxy2")
)

## With a model for Y_proxy_1 and Y_proxy_2
fit_lm_12 <- MLcovar(
  main_model = fn_true_lm,
  proxy_model = fn_proxy_lm,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE),
  args_proxy_model = list(method = "proxy_1_2")
)

## use_full = FALSE
fit_lm_1_ell <- MLcovar(
  main_model = fn_true_lm,
  proxy_model = fn_proxy_lm,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = FALSE),
  args_proxy_model = list(method = "proxy1")
)

## Estimation results
summary(fit_lm_1)
summary(fit_lm_2)
summary(fit_lm_12)
summary(fit_lm_1_ell)

## True value
1


## 2. Estimating the generalized linear regression coefficient
## Prepare an estimation model with true labels
fn_true_glm <- function (df) {
  fit <- glm(Y_binary ~ D + X1 + X2, family = binomial(link = "logit"), data = df)
  fit$coefficients["D"]
}

fn_proxy_glm <- function (df) {
  fit <- glm(Y_binary_proxy ~ D_proxy + X1_proxy + X2_proxy, family = binomial(link = "logit"), data = df)
  fit$coefficients
}

## Estimate with the proposed estimator
fit_glm <- MLcovar(
  main_model = fn_true_glm,
  proxy_model = fn_proxy_glm,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_glm)

## True value
2

## For the model with proxy variables, we can use other models than the model with labeled variables
fn_proxy_glm_lm <- function (df) {
  fit <- lm(Y_binary_proxy ~ D_proxy + X1_proxy + X2_proxy, data = df)
  fit$coefficients
}

fit_glm_lm <- MLcovar(
  main_model = fn_true_glm,
  proxy_model = fn_proxy_glm_lm,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_glm_lm)

## True value
2


## 3. Estimating the treatment effect by inverse probability weighting
## Prepare an estimation model with true labels
fn_true_ipw <- function (df) {
  fit_D <- glm(D ~ X1 + X2, family = binomial(link = "logit"), data = df)
  ps <- fit_D$fitted.values
  te <- sum(df$D * df$Y / ps) / sum(df$D / ps) - 
          sum((1 - df$D) * df$Y / (1 - ps)) / sum((1 - df$D) / (1 - ps))
  te
}

fn_proxy_ipw <- function (df) {
  fit_D <- glm(D_proxy ~ X1_proxy + X2_proxy, family = binomial(link = "logit"), data = df)
  ps <- fit_D$fitted.values
  te <- sum(df$D_proxy * df$Y_proxy_1 / ps) / sum(df$D_proxy / ps) - 
          sum((1 - df$D_proxy) * df$Y_proxy_1 / (1 - ps)) / sum((1 - df$D_proxy) / (1 - ps))
  te
}

## Estimate with the proposed estimator
fit_ipw <- MLcovar(
  main_model = fn_true_ipw,
  proxy_model = fn_proxy_ipw,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_ipw)

## True value
1


## 4. Estimating the treatment effect by augmented inverse probability weighting
## Prepare an estimation model with true labels
fn_true_aipw <- function (df) {
  fit_Y0 <- lm(Y ~ X1 + X2, data = subset(df, D == 0))
  pred_Y0 <- predict(fit_Y0, newdata = df)
  fit_Y1 <- lm(Y ~ X1 + X2, data = subset(df, D == 1))
  pred_Y1 <- predict(fit_Y1, newdata = df)
  fit_D <- glm(D ~ X1 + X2, family = binomial(link = "logit"), data = df)
  ps <- fit_D$fitted.values
  te <- mean(pred_Y1) - mean(pred_Y0) + 
          sum(df$D * (df$Y - pred_Y1) / ps) / sum(df$D / ps) - 
          sum((1 - df$D) * (df$Y - pred_Y0) / (1 - ps)) / sum((1 - df$D) / (1 - ps))
  te
}

fn_proxy_aipw <- function (df) {
  fit_Y0 <- lm(Y_proxy_1 ~ X1_proxy + X2_proxy, data = subset(df, D_proxy == 0))
  pred_Y0 <- predict(fit_Y0, newdata = df)
  fit_Y1 <- lm(Y_proxy_1 ~ X1_proxy + X2_proxy, data = subset(df, D_proxy == 1))
  pred_Y1 <- predict(fit_Y1, newdata = df)
  fit_D <- glm(D_proxy ~ X1_proxy + X2_proxy, family = binomial(link = "logit"), data = df)
  ps <- fit_D$fitted.values
  te <- mean(pred_Y1) - mean(pred_Y0) + 
          sum(df$D_proxy * (df$Y_proxy_1 - pred_Y1) / ps) / sum(df$D_proxy / ps) - 
          sum((1 - df$D_proxy) * (df$Y_proxy_1 - pred_Y0) / (1 - ps)) / sum((1 - df$D_proxy) / (1 - ps))
  te
}

## Estimate with the proposed estimator
fit_aipw <- MLcovar(
  main_model = fn_true_aipw,
  proxy_model = fn_proxy_aipw,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_aipw)

## True value
1


## 5. Estimating the local ATE with an instrumental variable
## Prepare an estimation model with true labels
fn_true_iv <- function (df) {
  fit <- estimatr::iv_robust(Y ~ D + X1 + X2 | Z + X1 + X2, data = df)
  fit$coefficients["D"]
}

fn_proxy_iv <- function (df) {
  fit <- estimatr::iv_robust(Y_proxy_1 ~ D_proxy + X1_proxy + X2_proxy | Z_proxy + X1_proxy + X2_proxy, data = df)
  fit$coefficients
}

## Estimate with the proposed estimator
fit_iv <- MLcovar(
  main_model = fn_true_iv,
  proxy_model = fn_proxy_iv,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_iv)

## True value
1

## We can use linear regression for proxy variables
fn_proxy_iv_lm <- function (df) {
  fit <- lm(Y_proxy_1 ~ D_proxy + X1_proxy + X2_proxy, data = df)
  fit$coefficients
}

fit_iv_lm <- MLcovar(
  main_model = fn_true_iv,
  proxy_model = fn_proxy_iv_lm,
  data = df,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_iv_lm)

## True value
1



## 6. Estimating the hazarad ratio with the Cox proportional hazard model
## Generate data
## Weibull distribution
## Latent event times (Weibull distribution)
lambda <- 1 / exp(-5 + 0.1 * (D + X1 + X2) / 5)
event <- rweibull(n = n_all, shape = 5, scale = lambda)

## Censoring times
censor <- pmin(runif(n_all, 50, 300), 250)

## Time
time <- pmin(event, censor)

## Status
status <- as.numeric(event < censor)

## Add noise to the treatment to create a proxy treatment (misclassifiction prob.: 3%)
status_proxy <- rbinom(n = n_all, size = 1, prob = 0.03 + status * 0.94)

# Create a data frame
df_survival <- data.frame(
  time = time,
  status = status,
  status_proxy = status_proxy,
  D = D,
  D_proxy = D_proxy,
  X1 = X1,
  X2 = X2,
  X1_proxy = X1_proxy,
  X2_proxy = X2_proxy,
  is_labeled = c(rep(1, n_ell), rep(0, n_u))
) |> 
  dplyr::mutate(status = dplyr::if_else(is_labeled == 1, status, NA_real_),
                D = dplyr::if_else(is_labeled == 1, D, NA_real_),
                X1 = dplyr::if_else(is_labeled == 1, X1, NA_real_),
                X2 = dplyr::if_else(is_labeled == 1, X2, NA_real_))

## Estimation
## Prepare an estimation model with true labels
fn_true_cox <- function (df) {
  fit <- survival::coxph(survival::Surv(time, status) ~ D + X1 + X2, data = df)
  fit$coefficients["D"]
}

fn_proxy_cox <- function (df) {
  fit <- survival::coxph(survival::Surv(time, status_proxy) ~ D_proxy + X1_proxy + X2_proxy, data = df)
  fit$coefficients
}

## Estimate with the proposed estimator
fit_cox <- MLcovar(
  main_model = fn_true_cox,
  proxy_model = fn_proxy_cox,
  data = df_survival,
  labeled_set_var_name = "is_labeled",
  options = SetOptions(n_boot = 5000, use_full = TRUE)
)
summary(fit_cox)

## True value
0.1

