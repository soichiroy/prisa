
#' @title Summary for MLcovar results
#' #' @description Provides a summary of the results from the MLcovar function.
#' 
#' @param object An object of class "MLcovar" returned by the MLcovar function.
#' @param ... Additional arguments.
#' @return A list containing the summary of the results.
#' @export
summary.MLcovar <- function(object, ...) {
  if (!inherits(object, "MLcovar")) {
    stop("The object must be of class 'MLcovar'.")
  }
  
  estimates <- object$estimates
  class(estimates) <- c(class(estimates), "summary.MLcovar")
  return(estimates)
}

#' @title Print summary for MLcovar results
#' 
#' @description Prints the summary of the results from the MLcovar function.
#' 
#' @param x An object of class "summary.MLcovar" returned by the summary.MLcovar
#'  function.
#' @param ... Additional arguments.
#' @export
print.summary.MLcovar <- function(x, ...) {
  cat("Estimates of the parameters:\n")
  print(as.data.frame(x))
  invisible(x)
}