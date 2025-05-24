#' @title Summary for MLcovar results
#' #' @description Provides a summary of the results from the MLcovar function.
#'
#' @param object An object of class "MLcovar" returned by the MLcovar function.
#' @param ... Additional arguments.
#' @return A list containing the summary of the results.
#' @export
summary.peri <- function(object, ...) {
  if (!inherits(object, "peri")) {
    stop("The object must be of class 'peri'.")
  }

  estimates <- object[c("estimates", "data_list")]
  class(estimates) <- c(class(estimates), "summary.peri")
  return(estimates)
}

#' @title Print summary for peri results
#'
#' @description Prints the summary of the results from the peri function.
#'
#' @param x An object of class "summary.peri" returned by the summary.peri
#'  function.
#' @param ... Additional arguments.
#' @importFrom cli cli_h1 cli_h2 cli_h3
#' @export
print.summary.peri <- function(x, digits = 4, ...) {
  cli::cli_h1("Prediction-error Robust Inference (peri) Results")
  cli::cli_h2("Main Estimates")
  print(
    format(as.data.frame(x$estimates$main), digits = digits),
    row.names = FALSE
  )
  cat("\n")

  # Show additional information
  cli::cli_h2("Additional Information")
  .PrintEfficiencyGain(
    estimates = x$estimates,
    data_list = x$data_list
  )

  .PrintDataInfo(x$data_list)

  cli::cli_h3("Labeled Only Estimates")
  print(
    format(as.data.frame(x$estimates$labeled_only), digits = digits),
    row.names = FALSE
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
    "Contribution from unlabeled set: {round(elss - n_ell, 1)}",
    "Standard Error Reduction: {std_err_reduction}%"
  ))

}
