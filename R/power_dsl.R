#' Power Analysis for DSL Regression
#' @param labeled_size A vector indicating the number of labeled documents for which the function predicts standard errors.
#' @param dsl_out An output from function \code{dsl}. When this is supplied, the remaining arguments are overwritten by arguments specified in the output of \code{dsl}. When this is \code{NULL}, the function will use arguments specified below.
#' @param model A regression model \code{dsl} currently supports \code{lm} (linear regression), \code{logit} (logistic regression), and \code{felm} (fixed-effects regression).
#' @param formula A formula used in the specified regression model.
#' @param predicted_var A vector of column names in the data that correspond to variables that need to be predicted.
#' @param prediction A vector of column names in the data that correspond to predictions of \code{predicted_var}.
#' @param data A data frame. The class should be \code{data.frame}.
#' @param cluster A column name in the data that indicates the level at which cluster standard errors are calculated. Default is \code{NULL}.
#' @param labeled (Optional) A column name in the data that indicates which observation is labeled. It should be a vector of 1 (labeled) and 0 (non-labeled). When \code{NULL}, the function assumes that observations that have \code{NA} in \code{predicted_var} are non-labeled and other observations are labeled.
#' @param sample_prob (Optional) A column name in the data that correspond to the sampling probability for labeling a particular observation. When \code{NULL}, the function assumes random sampling with equal probabilities.
#' @param fixed_effect (Used when \code{model = "felm"}) A type of fixed effects regression you run. \code{oneway} (one-way fixed effects) or \code{twoways} (two-way fixed effects).
#' @param index (Used when \code{model = "felm"}) A vector of column names specifying fixed effects. When \code{fixed_effect = oneway}, it has one element. When \code{fixed_effect = twoways}, it has two elements, e.g., \code{index = c("state", "year")}.
#' @param sl_method A name of a supervised machine learning model used internally to predict \code{predicted_var} by fine-tuning \code{prediction} or using predictors (specified in \code{feature}) when \code{prediction = NULL}. Users can run \code{available_method()} to see available supervised machine learning methods. Default is \code{grf} (generalized random forest).
#' @param feature A vector of column names in the data that correspond to predictors used to fit a supervised machine learning (specified in \code{sl_method}).
#' @param family (Used when making predictions) A variable type of \code{predicted_var}. Default is \code{gaussian}.
#' @param cross_fit The fold of cross-fitting. Default is \code{5}.
#' @param sample_split The number of sampling-splitting. Default is \code{10}.
#' @param seed Numeric \code{seed} used internally. Default is \code{1234}.
#' @return \code{dsl} returns an object of \code{dsl} class.
#'  \itemize{
#'    \item \code{predicted_se}: Predicted standard errors for coefficients. The first row shows the current standard errors for coefficients. The remaining rows show predicted standard errors.
#'    \item \code{labeled_size}: A vector indicating the number of labeled documents for which the function predicts standard errors.
#'    \item \code{dsl_out}: An output from function \code{dsl}.
#'  }
#' @export

power_dsl <- function(labeled_size = NULL,
                      dsl_out = NULL,
                      model = "lm",
                      formula,
                      predicted_var,
                      prediction = NULL,
                      data,
                      cluster = NULL,
                      labeled = NULL,
                      sample_prob = NULL,
                      index = NULL,
                      fixed_effect = "oneway",
                      sl_method = "grf",
                      feature = NULL,
                      family = "gaussian",
                      cross_fit = 5,
                      sample_split = 10,
                      seed = 1234){

  if(is.null(dsl_out) == TRUE){
    # Initial
    dsl_out <- dsl(model = model,
                   formula = formula,
                   predicted_var = predicted_var,
                   prediction = prediction,
                   labeled = labeled,
                   sample_prob = sample_prob,
                   feature = feature,
                   index = index,
                   fixed_effect = fixed_effect,
                   cluster = cluster,
                   data = data,
                   sl_method = sl_method,
                   family = family,
                   cross_fit = cross_fit,
                   sample_split = sample_split,
                   seed = seed)
  }

  sample_split_use <- dsl_out$internal$sample_split
  model_use <- dsl_out$internal$model
  full_se <- dsl_out$standard_errors
  ncol_X <- length(full_se)

  # Proportion
  full_size <- dsl_out$internal$num_expert
  if(is.null(labeled_size) == TRUE){
    labeled_size <- seq(from = 1.1, to = 2, by = 0.1)*full_size
  }
  pi_ratio <- (labeled_size/full_size)
  n <- dsl_out$internal$num_data

  cat("\nPower Analysis: ")
  exp_se0 <- array(NA, dim = c(length(pi_ratio), length(full_se), sample_split_use))
  for(ss_use in 1:sample_split_use){
    cat(paste0(ss_use,"/", sample_split_use, ".."))
    if(model_use == "felm"){
      s_J0  <- Matrix(dsl_out$internal$vcov_J[,,ss_use], sparse = TRUE)
      s_J <- solve(s_J0)
      s_J_use <- s_J[1:ncol_X, ,drop = FALSE]
    }else{
      s_J <- solve(dsl_out$internal$vcov_J[,,ss_use])
      s_J_use <- s_J[1:ncol_X, ,drop = FALSE]
    }

    for(z in 1:length(pi_ratio)){

      M0 <- (1/pi_ratio[z])*dsl_out$internal$vcov_main_1[,,ss_use] + dsl_out$internal$vcov_main_23[,,ss_use]

      if(model_use == "felm"){
        M0 <- Matrix(M0, sparse = TRUE)
      }

      V0 <- (s_J_use %*% M0 %*% t(s_J_use))/n
      V  <- dsl_out$internal$vcov_D[,,ss_use] %*% V0 %*% t(dsl_out$internal$vcov_D[,,ss_use])
      suppressWarnings({exp_se0[z, 1:ncol(V), ss_use] <- sqrt(diag(V))})
    }
  }

  exp_se <- apply(exp_se0, MARGIN = 1:2, median, na.rm = TRUE)

  rownames(exp_se) <- labeled_size
  colnames(exp_se) <- names(full_se)

  predicted_se <- rbind(full_se, exp_se)
  rownames(predicted_se)[1] <- full_size


  if(any(is.na(predicted_se))){
    # For some coefficients, predicted standard errors are too small and predictions are below zero.
    # For such coefficients, we use simple approximation
    pi_ratio_use <- c(1, pi_ratio)
    for(j in 1:ncol(predicted_se)){
      pos_ind <- 1
      for(i in seq(from = 2, to = nrow(predicted_se))){
        if(is.na(predicted_se[i, j])){
          predicted_se[i,j] <- predicted_se[pos_ind, j]/sqrt(pi_ratio_use[i]/pi_ratio_use[pos_ind])
        }else{
          pos_ind <- i
        }
      }
    }
  }

  out <- list("predicted_se" = predicted_se,
              "labeled_size" = c(full_size, labeled_size),
              "dsl_out" = dsl_out)

  class(out) <-c(class(out), "power_dsl")

  return(out)
}
