#' Estimating Regression using the DSL framework
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
#' @rawNamespace import(Matrix, except = summary)
#' @importFrom grf regression_forest
#' @importFrom estimatr lm_robust
#' @importFrom matrixcalc is.positive.definite
#' @importFrom arm model.matrixBayes
#' @importFrom stats as.formula glm lm median model.frame model.matrix model.response optim predict sd pnorm qnorm var
#' @importFrom utils capture.output
#' @importFrom graphics points
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
                sample_split = 10,
                seed = 1234){

  # ##################
  # Setup
  # ##################
  optim_method <-  "L-BFGS-B"
  lambda <-  0

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

  # labeled
  if(is.null(labeled) == TRUE){
    label_index <- 1 - as.numeric(apply(data[, predicted_var, drop = FALSE], 1, function(x) any(is.na(x))))
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
  covariates_use <- setdiff(covariates_use, predicted_var)

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
  if(any(is.na(data[data[, labeled] == 1, predicted_var]))){
    stop("NA in `predicted_var` of labeled data")
  }
  keep_var <- setdiff(unique(c("cluster___", sample_prob, labeled, covariates_use, all.vars(formula))), predicted_var)
  remove <- apply(data[, keep_var, drop = FALSE], 1, function(x) any(is.na(x)))
  data <- data[remove == FALSE, , drop = FALSE]

  # compute the total num of samples and labeled data
  num_expert <- sum(data[, labeled], na.rm = TRUE)
  num_data   <- nrow(data)

  # Standardize sample_prob
  mean_prob <- sum(data[, labeled], na.rm = TRUE)/nrow(data)
  data[, sample_prob] <- data[, sample_prob]*(mean_prob/mean(data[, sample_prob], na.rm = TRUE))

  if(any(is.na(data[data[, labeled] == 0, predicted_var]) == FALSE)){
    stop("Some `predicted_var` in non-labeled data are not NA. Please check the data. ")
  }

  # Checking values
  if(all(sort(unique(data[, labeled])) %in% c(0,1)) == FALSE){
    stop(" `labeled` in `data` should be 0 or 1. Please check the data. ")
  }
  if(all(data[, sample_prob] > 0 & data[, sample_prob] <= 1) == FALSE){
    stop(" `sample_prob` in `data` should be greater than 0 and equal to or smaller than 1. Please check the data. ")
  }


  # ##################
  # Data Preparation
  # ##################
  # Unique Cluster
  uniq_cluster <- sort(unique(data$cluster___))

  # Number of Columns we analyze
  ncol_X <- ncol(model.matrix(formula, data = data))
  if(model == "felm"){
    ncol_X <- ncol_X - 1   # Remove Intercept

    ## including fixed effects
    if(fixed_effect != "twoways"){
      ncol_X_exp <- ncol_X + length(unique(data[, index_use]))
    }else if(fixed_effect == "twoways"){
      ncol_X_exp <- ncol_X + length(unique(data[, index[1]])) + length(unique(data[, index[2]])) - 1
    }
  }else{
    ncol_X_exp <- ncol_X
  }

  # matrices
  dm_point_l <- dm_se_l <- matrix(NA, ncol = ncol_X, nrow = sample_split)
  dm_vcov_l <- array(NA, dim = c(ncol_X, ncol_X, sample_split))
  dm_vcov_main_1_l <- dm_vcov_main_23_l <- array(NA, dim = c(ncol_X_exp, ncol_X_exp, sample_split))
  dm_vcov_J_l <- array(NA, dim = c(ncol_X_exp, ncol_X_exp, sample_split))
  dm_vcov_D_l <- array(NA, dim = c(ncol_X, ncol_X, sample_split))
  RMSE_cv <- matrix(NA, nrow = length(predicted_var), ncol = sample_split)

  # ######################
  # Sample Splitting
  # ######################

  set.seed(seed)

  cat("Sample Splitting: ")
  for(ss_use in 1:sample_split){
    cat(paste0(ss_use,"/", sample_split, ".."))

    id_base0 <- rep(1:cross_fit, each = floor(length(uniq_cluster)/cross_fit))
    extra <- length(uniq_cluster) - length(id_base0)
    if(extra > 0){
      id_base0 <- c(id_base0, 1:extra)
    }else{
      id_base  <- id_base0
    }
    id_base  <- sample(id_base0, size = length(id_base0), replace = FALSE)

    dm_predicted_var <- list()
    id_final <- c()
    RMSE_cv0 <- matrix(NA, nrow = length(predicted_var), ncol = cross_fit)
    # ############
    # Cross-Fit
    # ############
    for(i_use in 1:cross_fit){
      train_id <- uniq_cluster[id_base != i_use]
      test_id  <- setdiff(uniq_cluster, train_id)

      train_which <- unlist(sapply(train_id, function(x) which(data$cluster___ == x)))
      test_which  <- unlist(sapply(test_id, function(x) which(data$cluster___ == x)))

      data_train     <- data[train_which, , drop = FALSE]
      data_test      <- data[test_which, , drop = FALSE]

      # Fit the model E( predicted_var | covariates_use, R = 1) in the training sample

      for(v_use in 1:length(predicted_var)){
        fit_train <- fit_model(outcome = predicted_var[v_use], labeled = labeled, prediction = prediction,
                               covariates = covariates_use,
                               data = data_train[data_train[, labeled] == 1, , drop = FALSE],
                               seed = seed, method = sl_method, family = family)

        # Predict E(predicted_var | covariates_use, R = 1) in the test data
        out_test <- fit_test(fit_out = fit_train,
                             outcome = predicted_var[v_use], labeled = labeled, covariates = covariates_use,
                             method = sl_method, family = family,
                             data = data_test, seed = seed)

        predicted_var_test <- out_test$Y_hat

        # Store results
        if(i_use == 1){
          dm_predicted_var[[v_use]] <- predicted_var_test
        }else{
          dm_predicted_var[[v_use]] <- c(predicted_var_test, dm_predicted_var[[v_use]])
        }

        RMSE_cv0[v_use, i_use] <- out_test$RMSE
      }
      rm(data_train)
      rm(data_test)
      # id
      id_final <- c(test_which, id_final)
    }
    rownames(RMSE_cv0) <- predicted_var

    # data we use for estimation
    data_pred <- data_orig <- data[id_final, , drop = FALSE]

    data_pred[, predicted_var] <- do.call("cbind", dm_predicted_var)

    # ##############
    # DSL Estimator
    # ##############
    fit_dm <- dsl_general_moment_est(model = model,
                                     formula = formula,
                                     labeled = labeled,
                                     sample_prob = sample_prob,
                                     predicted_var = predicted_var,
                                     data_orig = data_orig,
                                     data_pred = data_pred,
                                     index = index,
                                     fixed_effect = fixed_effect,
                                     clustered = clustered,
                                     optim_method = optim_method,
                                     lambda = lambda)

    dm_point_l[ss_use, 1:ncol_X] <- fit_dm$coefficients
    dm_se_l[ss_use, 1:ncol_X]    <- fit_dm$standard_errors
    dm_vcov_l[1:ncol_X, 1:ncol_X, ss_use] <- as.matrix(fit_dm$vcov)

    # dm_vcov_meat_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$Meat)
    dm_vcov_main_1_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$Meat_decomp$main_1)
    dm_vcov_main_23_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$Meat_decomp$main_23)
    dm_vcov_J_l[1:ncol_X_exp, 1:ncol_X_exp, ss_use] <- as.matrix(fit_dm$J)
    dm_vcov_D_l[1:ncol_X, 1:ncol_X, ss_use] <- fit_dm$D

    RMSE_cv[1:length(predicted_var), ss_use] <- apply(RMSE_cv0, 1, mean, na.rm = TRUE)
  }
  rm(data)

  # ################################
  # Combine across sampling-splitting (Take into account the variation across sample splitting)
  # ################################

  # point estimates + standard errors
  dm_point <- apply(dm_point_l, 2, median)
  dm_se_0  <- dm_se_l^2 + (t(t(dm_point_l) - dm_point))^2
  dm_se    <- sqrt(apply(dm_se_0, 2, median))
  names(dm_point) <- names(dm_se) <- names(fit_dm$coef)

  # vcov (Take into account the variation across sample splitting)
  dm_vcov_mat <- array(NA, dim = c(ncol_X, ncol_X, sample_split))
  for(ss_use in 1:sample_split){
    dm_point_use <- dm_point_l[ss_use, 1:ncol_X] - dm_point
    dm_vcov_adj  <- matrix(dm_point_use, nrow = ncol_X, ncol =1) %*% matrix(dm_point_use, nrow = 1, ncol = ncol_X)
    dm_vcov_mat[1:ncol_X, 1:ncol_X, ss_use] <- dm_vcov_l[1:ncol_X, 1:ncol_X, ss_use] + dm_vcov_adj
  }
  dm_vcov0 <- apply(dm_vcov_mat, MARGIN = 1:2, median)
  dm_vcov  <- (dm_vcov0 + t(dm_vcov0))/2 # make it symmetric
  if(is.positive.definite(dm_vcov) == FALSE){
    dm_vcov <- nearPD(dm_vcov, keepDiag = TRUE)
  }

  internal <- list("formula" = formula, "model" = model,
                   "predicted_var" = predicted_var, "prediction" = prediction,
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
              "vcov" = dm_vcov,
              "RMSE" = RMSE_cv,
              "internal" = internal)

  class(out)  <- c(class(out), "dsl")
  return(out)
}
