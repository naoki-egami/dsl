#' Estimating Regression using the DSL framework
#' @param model A regression model \code{dsl} currently supports \code{lm} (linear regression), \code{logit} (logistic regression), and \code{felm} (fixed-effects regression).
#' @param formula A formula used in the specified regression model.
#' @param predicted_var A vector of column names in the data that correspond to variables that need to be predicted. When \code{unc_label = TRUE}, this should be \code{list} where each element of the list contains variable names specifying labels from different coders.
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
#' @param unc_label A logical value (\code{True}, or \code{FALSE}) indicating whether to apply the quasi-Bayesian approach to incorporate uncertainties in labels. Default is \code{FALSE}.
#' @param unc_label_sim (Used when \code{unc_label = TRUE}) The number of re-sampling of labels within each observation. Default is \code{100}.
#' @param tuning A logical value (\code{True}, or \code{FALSE}) indicating whether to apply the power-tuning. Default is \code{TRUE}.
#' @param seed Numeric \code{seed} used internally. Default is \code{1234}.
#' @rawNamespace import(Matrix, except = summary)
#' @importFrom grf regression_forest
#' @importFrom estimatr lm_robust
#' @importFrom matrixcalc is.positive.definite
#' @importFrom arm model.matrixBayes
#' @importFrom stats as.formula glm lm median model.frame model.matrix model.response optim predict sd pnorm qnorm var cov
#' @importFrom utils capture.output
#' @importFrom graphics points
#' @importFrom MASS mvrnorm
#' @import tidyverse
#' @import SuperLearner
#' @return \code{dsl} returns an object of \code{dsl} class.
#'  \itemize{
#'    \item \code{coefficients}: Estimated coefficients.
#'    \item \code{standard_errors}: Estimated standard errors.
#'    \item \code{vcov}: Estimated variance-covariance matrix.
#'    \item \code{RMSE}: Root mean squared error in the internal prediction step.
#'    \item \code{internal}: Outputs used only for the internal use.
#'  }
#' @export

dsl <- function(model = "lm",
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
                sample_split = 5,
                unc_label = FALSE,
                unc_label_sim = 50,
                tuning = TRUE,
                seed = 1234){

  # ##################
  # Setup
  # ##################
  optim_method <-  "L-BFGS-B"
  lambda <-  0
  tuning_para <- 1

  # data.frame
  class(data) <- "data.frame"

  if((model %in% c("lm", "logit", "felm")) == FALSE){
    stop(" `model` should be either `lm`, `logit`, or `felm` ")
  }

  if(is.null(prediction) & is.null(feature)){
    stop(" Please specify at least one of `prediction` and `feature` ")
  }

  # Model-Specific Check
  if(model == "felm"){
    if(is.null(index) == TRUE){
      stop("Please specify `index` for `model = felm` ")
    }
    if(fixed_effect == "twoways"){
      if(length(index) != 2){
        stop("Please specify two indices for the two-way fixed effects estimator.")
      }
    }else if(fixed_effect == "oneway"){
      index_use <- index[1]
    }
  }

  # predicted_var
  if(unc_label == TRUE){
    if(is.list(predicted_var) == FALSE){
      stop("`predicted_var` should be a `list` when `unc_label = TRUE`.")
    }
    predicted_var_all <- unlist(predicted_var)
  }else{
    predicted_var_all <- predicted_var
  }

  # labeled
  if(is.null(labeled) == TRUE){
    label_index <- 1 - as.numeric(apply(data[, predicted_var_all, drop = FALSE], 1, function(x) any(is.na(x))))
    data$labeled___ <- label_index
    labeled <- "labeled___"

    if(mean(data[,labeled]) == 1){
      warnings(" Please use `NA` to indicate observations that do not have labels. ")
    }
  }

  # sample_prob
  if(is.null(sample_prob) == TRUE){
    sample_prob0 <- sum(data[, labeled], na.rm = TRUE)/nrow(data)
    data$sample_prob___ <- sample_prob0
    sample_prob <- "sample_prob___"
    equal_prob <- TRUE
  }else{
    if(var(data[, sample_prob]) == 0){
      equal_prob <- TRUE
    }else{
      equal_prob <- FALSE
    }
  }

  # Setup feature
  covariates_use <- unique(c(prediction, feature))
  covariates_use <- setdiff(covariates_use, predicted_var_all)

  # Handle Methods
  if(any(sl_method != "grf")){
    all_method <- available_method(print_out = FALSE)
    if(any((sl_method %in% all_method) == FALSE)){
      stop(" `sl_method` is not supported in `dsl`. Please run `available_method() to see what `sl_method` is supported. ")
    }
  }

  # Handle cluster standard errors
  if(is.null(cluster) == TRUE){
    data$cluster___ <- seq(1:nrow(data))
    clustered <- FALSE
  }else{
    data$cluster___ <- as.numeric(as.factor(data[, cluster]))
    clustered <- TRUE
  }

  # Handling NAs
  if(any(is.na(data[data[, labeled] == 1, predicted_var_all]))){
    stop("NA in `predicted_var` of labeled data")
  }
  all_var_for <- all.vars(formula)
  if(unc_label == TRUE){
    all_var_for <- setdiff(all_var_for, names(predicted_var))
  }
  keep_var <- setdiff(unique(c("cluster___", sample_prob, labeled, covariates_use, all_var_for )), predicted_var_all)
  remove <- apply(data[, keep_var, drop = FALSE], 1, function(x) any(is.na(x)))
  data <- data[remove == FALSE, , drop = FALSE]

  # compute the total num of samples and labeled data
  num_expert <- sum(data[, labeled], na.rm = TRUE)
  num_data   <- nrow(data)

  # Standardize sample_prob
  mean_prob <- sum(data[, labeled], na.rm = TRUE)/nrow(data)
  data[, sample_prob] <- data[, sample_prob]*(mean_prob/mean(data[, sample_prob], na.rm = TRUE))

  if(any(is.na(data[data[, labeled] == 0, predicted_var_all]) == FALSE)){
    stop("Some `predicted_var` in non-labeled data are not NA. Please check the data. ")
  }

  # Checking values
  if(all(sort(unique(data[, labeled])) %in% c(0,1)) == FALSE){
    stop(" `labeled` in `data` should be 0 or 1. Please check the data. ")
  }
  if(all(data[, sample_prob] > 0 & data[, sample_prob] <= 1) == FALSE){
    stop(" `sample_prob` in `data` should be greater than 0 and equal to or smaller than 1. Please check the data. ")
  }


  # #####################################
  # Uncertainties in Expert-Annotations (unc_label = TRUE)
  # #####################################
  if(unc_label == FALSE){
    predicted_var_use <- predicted_var

    # if(tuning == TRUE){
    #   est_tuning_para <- fit_dsl(data, formula, model, fixed_effect, index_use, index,
    #                              predicted_var = predicted_var_use, labeled, prediction, covariates_use,
    #                              sample_prob, clustered, feature,
    #                              seed, sl_method, family, sample_split = 1, cross_fit,
    #                              num_expert, equal_prob, num_data, cluster,
    #                              optim_method, lambda, verbose = FALSE, tuning = TRUE, tuning_para = 1)
    # }else{
    #   est_tuning_para <- 1
    # }

    out <- fit_dsl(data, formula, model, fixed_effect, index_use, index,
                   predicted_var = predicted_var_use, labeled, prediction, covariates_use,
                   sample_prob, clustered, feature,
                   seed, sl_method, family, sample_split, cross_fit,
                   num_expert, equal_prob, num_data, cluster,
                   optim_method, lambda, verbose = TRUE, tuning = tuning, tuning_para = tuning_para)

  }else if(unc_label == TRUE){
    # randomly sampling
    predicted_var_use <- names(predicted_var)
    predicted_var_data <- data.frame(matrix(NA, nrow = nrow(data), ncol = length(predicted_var_use)))
    colnames(predicted_var_data) <- predicted_var_use
    RMSE_cv_M <- matrix(NA, nrow = length(predicted_var_use), ncol = unc_label_sim)

    data_use <- data[, is.element(colnames(data), predicted_var_use) == FALSE]

    set.seed(seed)
    cat("Re-Sampling Labels: ")
    for(rs in 1:unc_label_sim){
      perct <- round(rs*100/unc_label_sim)
      if(perct%%10 == 0){
        cat(paste0(perct, "%.."))
      }

      # re-sampling
      for(num_pr in 1:length(predicted_var_use)){
        use_exp <- sample(x = seq(1:length(predicted_var[[num_pr]])), size = nrow(data), replace = TRUE)
        for(num_exp in 1:length(predicted_var[[num_pr]])){
          predicted_var_data[use_exp == num_exp, num_pr] <- data[use_exp == num_exp, predicted_var[[num_pr]][num_exp]]
        }
      }
      data_rs <- cbind(data_use, predicted_var_data)

      fit_rs <- fit_dsl(data = data_rs, formula, model, fixed_effect, index_use, index,
                        predicted_var = predicted_var_use, labeled, prediction, covariates_use,
                        sample_prob, clustered, feature,
                        seed, sl_method, family, sample_split, cross_fit,
                        num_expert, equal_prob, num_data, cluster,
                        optim_method, lambda, verbose = FALSE, tuning = tuning)

      coef_rs_b <- mvrnorm(n = 100, mu = fit_rs$coefficients, Sigma = fit_rs$vcov) # nrow = 100, ncol = length(coefficient)
      RMSE_cv_M[1:nrow(RMSE_cv_M), rs]  <- apply(fit_rs$RMSE, 1, mean)

      if(rs == 1){
        coef_rs  <- coef_rs_b
        se_check <- fit_rs$standard_errors
      }else{
        coef_rs <- rbind(coef_rs, coef_rs_b)
        se_check <- rbind(se_check, fit_rs$standard_errors)
      }
    }

    # ##########
    # summarize
    # ##########
    dm_point <- apply(coef_rs, 2, mean, na.rm = TRUE)
    dm_se    <- apply(coef_rs, 2, sd, na.rm = TRUE)
    dm_vcov   <- cov(coef_rs, use = "complete.obs")
    RMSE_cv  <- apply(RMSE_cv_M, 1, mean, na.rm = TRUE)

    # for now
    dm_vcov_main_1_l <- NULL
    dm_vcov_main_23_l <- NULL
    dm_vcov_J_l <- NULL
    dm_vcov_D_l <- NULL

    internal <- list("formula" = formula, "model" = model,
                     "predicted_var" = predicted_var_use, "prediction" = prediction,
                     "num_expert" = num_expert, "equal_prob" = equal_prob, "num_data" = num_data,
                     "cluster" = cluster,
                     "index" = index, "fixed_effect" = fixed_effect,
                     "sl_method"= sl_method, "feature" = feature, "cross_fit" = cross_fit, "sample_split" = sample_split,
                     "vcov_main_1" = dm_vcov_main_1_l,
                     "vcov_main_23" = dm_vcov_main_23_l,
                     "vcov_J" = dm_vcov_J_l,
                     "vcov_D" = dm_vcov_D_l)

    out <- list("coefficients" = dm_point,
                "standard_errors" = dm_se,
                "se_check" = se_check,
                "vcov" = dm_vcov,
                "RMSE" = RMSE_cv,
                "internal" = internal)
  }

  class(out)  <- c(class(out), "dsl")
  return(out)
}
