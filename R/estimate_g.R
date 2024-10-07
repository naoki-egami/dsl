# Fit model
fit_model <- function(outcome, labeled, covariates, data, method, sample_prob = NULL, prediction = NULL, family = "gaussian", seed = 1234){

  formula   <- as.formula(paste0(outcome, "~", paste(c(labeled, covariates), collapse = "+")))
  formula_X <- as.formula(paste0(outcome, "~", paste(c(covariates), collapse = "+")))

  mf <- model.frame(formula, data = data)
  Y  <- as.numeric(model.response(mf))
  R <- as.numeric(mf[, labeled])
  X_name <- colnames(model.matrix(formula_X, data = data)) # this is more robust way to handle categorical variables
  X  <- model.matrix(formula, data = data)[, X_name[-1], drop = FALSE]   # this is more robust way to handle categorical variables

  set.seed(seed)
  if(any(method == "grf")){
    if(is.null(sample_prob) == FALSE){
      fit_out <- regression_forest(X = X[R == 1, , drop = FALSE], Y = Y[R == 1], seed = seed, sample.weights = 1/data[,sample_prob])
    }else if(is.null(sample_prob) == TRUE){
      fit_out <- regression_forest(X = X[R == 1, , drop = FALSE], Y = Y[R == 1], seed = seed)
    }
  }else if(method == "ranger"){
    method <- setdiff(method, "grf")
    sl_method <- paste0("SL.", method)
    X <- data.frame(X)
    suppressWarnings({fit_out <- SuperLearner::SuperLearner(Y = Y[R == 1], X = X[R == 1, , drop = FALSE],
                                                            SL.library = sl_method, family = family,
                                                            always.split.variables = prediction)})
  }else if(any(method == "identity")){
    fit_out <- "identity_function"
  }else if(any(method != "grf")){
    method <- setdiff(method, "grf")
    sl_method <- paste0("SL.", method)
    X <- data.frame(X)
    suppressWarnings({fit_out <- SuperLearner::SuperLearner(Y = Y[R == 1], X = X[R == 1, , drop = FALSE],
                                                            SL.library = sl_method, family = family)})
  }
  return(fit_out)
}

# Predict based on fitted g model
fit_test <- function(fit_out, outcome, labeled, covariates, data, method, family, seed = 1234){

  data$internal_id_for_test <- seq(1:nrow(data))
  formula_X <- as.formula(paste0("~ ", paste(c("internal_id_for_test", covariates), collapse = "+")))
  X_all  <- model.matrix(formula_X, data = data)[, -1, drop = FALSE]
  X <- X_all[, c(-1), drop = FALSE]

  # outcomes
  formula   <- as.formula(paste0(outcome, "~", paste(c("internal_id_for_test", labeled, covariates), collapse = "+")))
  mf <- model.frame(formula, data = data)
  Y  <- as.numeric(model.response(mf))

  if(any(method == "grf")){
    new_data_use <- X
    Y_hat <- predict(fit_out, newdata = new_data_use)$predictions
  }else if(any(method == "identity")){
    Y_hat <- as.numeric(X[,1])
  }else if(any(method != "grf")){
    new_data_use <- as.data.frame(X)
    suppressWarnings({Y_hat <- predict(fit_out, newdata = new_data_use, onlySL = TRUE)$pred})
  }

  # Cross-validation error
  Y_hat_use <- Y_hat[match(mf$internal_id_for_test, X_all[, "internal_id_for_test"])]
  RMSE <- sqrt(mean((Y - Y_hat_use)^2, na.rm = TRUE))

  out <- list("Y_hat" = Y_hat, "RMSE" = RMSE)

  return(out)
}

#' Showing available methods for Supervised Machine Learning
#' @param print_out Whether to print the output (\code{TRUE}) or not (\code{FALSE}).
#' @importFrom SuperLearner listWrappers
#' @return \code{availabel_method} returns a list of available methods.
#' @export
#'
available_method <- function(print_out = TRUE){
  suppressMessages(capture.output(sl_method <- listWrappers(what = "SL"), file = nullfile()))
  sl_method <- sl_method[regexpr("SL.", sl_method) == 1]
  sl_method <- gsub("SL.", "", sl_method)
  sl_method <- setdiff(sl_method, "template")

  all_method <- c("grf", sl_method, "identity")
  if(print_out == TRUE){
    print(all_method)
  }
  invisible(all_method)
}
