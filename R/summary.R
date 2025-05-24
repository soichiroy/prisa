
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
  
  estimates <- object$estimates
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
#' @export
print.summary.peri <- function(x, ...) {
  cat("Estimates of the parameters:\n")
  print(as.data.frame(x))
  invisible(x)
}