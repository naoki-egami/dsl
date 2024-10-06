fit_dsl <- function(data, formula, model, fixed_effect, index_use, index,
                    predicted_var, labeled, prediction, covariates_use,
                    sample_prob, clustered, feature,
                    seed, sl_method, family, sample_split, cross_fit,
                    num_expert, equal_prob, num_data, cluster,
                    optim_method, lambda, verbose) {

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

  if(verbose == TRUE){cat("Cross-Fitting: ")}
  for(ss_use in 1:sample_split){
    if(verbose == TRUE){cat(paste0(ss_use,"/", sample_split, ".."))}

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

  return(out)
}
