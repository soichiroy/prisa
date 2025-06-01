#' @title Summary for prisa results
#' @description Provides a summary of the results from the prisa function.
#'
#' @param object An object of class "prisa" returned by the prisa function.
#' @param ... Additional arguments.
#' @return A list containing the summary of the results.
#' @export
summary.prisa <- function(object, ...) {
  if (!inherits(object, "prisa")) {
    stop("The object must be of class 'prisa'.")
  }

  estimates <- object[c("estimates", "data_list")]
  class(estimates) <- c(class(estimates), "summary.prisa")
  return(estimates)
}

#' @title Print summary for prisa results
#'
#' @description Prints the summary of the results from the prisa function.
#'
#' @param x An object of class "summary.prisa" returned by the summary.prisa
#'  function.
#' @param digits The number of significant digits to print. Default is 4.
#' @param ... Additional arguments.
#' @importFrom cli cli_h1 cli_h2 cli_h3
#' @export
print.summary.prisa <- function(x, digits = 4, ...) {
  cli::cli_h1("Prediction-error Robust Inference (prisa) Results")
  cli::cli_h2("Main Estimates")
  print(
    format(as.data.frame(x$estimates$main), digits = digits)
  )
  cat("\n\n")
  # Show additional information
  cli::cli_h2("Additional Information")
  .PrintEfficiencyGain(
    estimates = x$estimates,
    data_list = x$data_list
  )

  .PrintDataInfo(x$data_list)

  cli::cli_h3("Labeled Only Estimates")
  print(
    format(as.data.frame(x$estimates$labeled_only), digits = digits)
  )
  invisible(x)
}

#' @title Prepare data information for printing
#' @importFrom cli cli_h3 cli_ul
#' @noRd
.PrintDataInfo <- function(data_list) {
  n_ell <- data_list$n_ell
  n_full <- data_list$n_full
  prop <- round(data_list$prop, 4)

  cli::cli_h3("Data")
  cli::cli_ul(c(
    "n_labeled_data: {n_ell}",
    "n_full_data: {n_full}",
    "proportion_of_coded_obs: {prop}"
  ))
  invisible(data_list)
}

#' @title Print efficiency gain from PERI
#' @importFrom cli cli_h3 cli_ul
#' @noRd
.PrintEfficiencyGain <- function(estimates, data_list) {
  std_err_reduction <- estimates$main$std_err / estimates$labeled_only$std_err
  std_err_reduction <- round(100 * (std_err_reduction - 1))
  elss <- estimates$main$elss
  n_ell <- data_list$n_ell

  cli::cli_h3("Efficiency Gain from PERI")
  cli::cli_ul(c(
    "Effective Labeled Sample Size (ELSS): {round(elss, 3)}",
    "Contribution from the unlabeled set: {round(elss - n_ell, 1)}",
    "Standard Error Reduction: {std_err_reduction}%"
  ))
  invisible(estimates)
}

#' @title Obtain estimation results
#'
#' @param x An object of class "peri" returned by the peri function.
#' @return A tibble with the following columns:
#' \describe{
#'   \item{estimator}{The type of estimator used (Prediction-error robust
#'    inference (peri) or labeled_only).}
#'   \item{estimate}{The estimated value.}
#'   \item{std_err}{The standard error of the estimate.}
#'   \item{ci_lower_95}{The lower bound of the 95% confidence interval.}
#'   \item{ci_upper_95}{The upper bound of the 95% confidence interval.}
#'   \item{elss}{The effective labeled sample size. For the labeled_only
#'    estimator, this corresponds to the number of labeled observations. For the
#'    peri, the estimated ELSS is shown.}
#' }
#' @seealso [peri()]
#' @export
#' @importFrom dplyr mutate bind_rows select everything
#' @importFrom tibble as_tibble rownames_to_column
#' @importFrom rlang .data
get_estimates <- function(x) {
  if (!inherits(x, "peri")) {
    stop("The object must be of class 'peri'.")
  }

  df_main <- rownames_to_column(x$estimates$main, var = "variable")
  df_labeled_only <- rownames_to_column(
    x$estimates$labeled_only,
    var = "variable"
  )
  out <- bind_rows(
    mutate(df_main, estimator = "peri"),
    mutate(df_labeled_only, estimator = "labeled_only")
  ) %>%
    select(.data$variable, .data$estimator, everything())
  as_tibble(out)
}
