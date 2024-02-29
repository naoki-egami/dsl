#' Summary function for \code{dsl}
#' @param object An output from function \code{dsl}.
#' @param ci A coverage rate of confidence intervals. Default is \code{0.95}, which reports 95\% confidence intervals.
#' @param digits The number of digits reported.
#' @param ... Other arguments.
#' @export

summary.dsl <- function(object, ci = 0.95, digits = 4, ...){

  out_tab <- cbind(object[[1]], object[[2]])
  alpha_h <- 1 - (1 - ci)/2
  ci_low  <- object[[1]] - qnorm(alpha_h)*object[[2]]
  ci_high <- object[[1]] + qnorm(alpha_h)*object[[2]]
  p_v <- 1 - pnorm(abs(object[[1]]/object[[2]]))
  out_tab <- cbind(object[[1]], object[[2]], ci_low, ci_high, p_v)
  colnames(out_tab) <- c("Estimate", "Std. Error", "CI Lower", "CI Upper", "p value")
  if(nrow(out_tab) == 1){
    rownames(out_tab) <- ""
  }else{
    rownames(out_tab) <- names(object[[1]])
  }

  # add significance
  sig <- rep("", length(out_tab[,5]))
  sig[out_tab[,5]  < 0.001] <- "***"
  sig[out_tab[,5]  >= 0.001 & out_tab[,5]  < 0.01] <- "**"
  sig[out_tab[,5]  >= 0.01 & out_tab[,5]  < 0.05] <- "*"
  sig[out_tab[,5]  >= 0.05 & out_tab[,5]  < 0.1] <- "."
  out_tab_print <- as.data.frame(out_tab)
  out_tab_print$Sig <- sig
  colnames(out_tab_print)[6] <- ""

  out_tab_print[, 1:5] <- round(out_tab_print[, 1:5], digits = digits)

  ch_for <- as.character(object$internal$formula)

  cat("==================\n")
  cat("DSL Specification:\n")
  cat("==================\n")

  if(object$internal$model != "felm"){
    cat(paste0("Model:  ", object$internal$model))
    cat("\n")
  }else if(object$internal$model == "felm"){
    cat(paste0("Model:  ", object$internal$model, " (", object$internal$fixed_effect, ")"))
    cat("\n")
  }

  cat(paste0("Call:  ", paste(ch_for[2], ch_for[1], ch_for[3])))
  cat("\n")

  if(object$internal$model == "felm"){
    if(object$internal$fixed_effect == "twoways"){
      cat(paste0("Fixed Effects:  ", object$internal$index[1], " and ", object$internal$index[2]))
      cat("\n\n")
    }else if(object$internal$fixed_effect == "oneway"){
      cat(paste0("Fixed Effects:  ", object$internal$index[1]))
      cat("\n\n")
    }
  }else{
    cat("\n")
  }

  cat(paste0("Predicted Variables:  ", paste(c(object$internal$predicted_var), collapse = ", ")))
  cat("\n")

  if(is.null(object$internal$prediction) == FALSE){
    cat(paste0("Prediction:  ", paste(c(object$internal$prediction), collapse = ", ")))
    cat("\n\n")
  }else{
    cat(paste0("Prediction: ", object$internal$method, " with predictors specified in `covariates`"))
    cat("\n\n")
  }

  cat(paste0("Number of Labeled Observations:  ", object$internal$num_expert))
  cat("\n")

  if(object$internal$equal_prob == TRUE){
    cat("Random Sampling for Labeling with Equal Probability: Yes")
    cat("\n\n")
  }else{
    cat("Random Sampling for Labeling with Equal Probability:  No\n")
    cat("(Sampling probabilities are defined in `sample_prob`)")
    cat("\n\n")
  }

  cat("=============\n")
  cat("Coefficients:\n")
  cat("=============\n")
  print(out_tab_print, row.names = TRUE)
  cat("---\n")
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  cat(paste0("\n", round(ci*100), "% confidence intervals (CI) are reported."))
  if(is.null(object$internal$cluster) == FALSE){
    cat(paste0("\nStandard errors are clustered by ", object$internal$cluster, "."))
  }

  invisible(out_tab_print)
}
