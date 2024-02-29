#' Summarizing estimates from DSL estimators
#' @param object An output from function \code{power_dsl}.
#' @param ... Other arguments.
#' @export

summary.power_dsl <- function(object,...){

  cat("\n")
  cat("==========================\n")
  cat("Predicted Standard Errors:\n")
  cat("==========================\n")
  print(object$predicted_se, row.names = TRUE)
  cat("---\n")

  invisible(object$predicted_se)
}
