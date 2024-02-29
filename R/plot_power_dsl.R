#' Summarizing estimates from DSL estimators
#' @param x An output from function \code{power_dsl}.
#' @param coef_name A name of a coefficient users want to visualize. When \code{NULL}, the function visualizes every coefficient.
#' @param ... Other arguments.
#' @importFrom graphics par
#' @export

plot.power_dsl <- function(x, coef_name = NULL, ...){

  if(is.null(coef_name) == TRUE){
    coef_name <- colnames(x$predicted_se)
  }

  if(length(coef_name) > 1){
    op <- par(ask=TRUE)
    for(z in 1:length(coef_name)){
      plot(x$labeled_size, x$predicted_se[,z], pch = 19, type = "b",
           xlab = "Size of Labeled Data", ylab = "Standard Error",
           main = coef_name[z], col = "red")
      points(x$labeled_size[1], x$predicted_se[1,z], pch = 15, cex = 1.5, col = "black")
    }
    par(op)
  }else if(length(coef_name) == 1){
    pick <- which(colnames(x$predicted_se) == coef_name)
    plot(x$labeled_size, x$predicted_se[,pick], pch = 19, type = "b",
         xlab = "Size of Labeled Data", ylab = "Standard Error",
         main = coef_name, col = "red")
    points(x$labeled_size[1], x$predicted_se[1,pick], pch = 15, cex = 1.5, col = "black")
  }
}
