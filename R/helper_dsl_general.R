dsl_general_moment_est <- function(model, formula, labeled, sample_prob, predicted_var, data_orig, data_pred, index = NULL, fixed_effect = NULL,
                                   clustered, optim_method = "L-BFGS-B",
                                   lambda = 0.00001, tuning, tuning_para){

  # Prepare data
  labeled_ind <- data_orig[, labeled]
  sample_prob_use <- data_orig[, sample_prob]

  mf_orig <- model.frame(formula, data = data_orig, na.action = "na.pass")
  X_orig <- model.matrix(formula, data = mf_orig)
  Y_orig <- as.numeric(model.response(mf_orig))
  mf_pred <- model.frame(formula, data = data_pred, na.action = "na.pass")
  X_pred <- model.matrix(formula, data = mf_pred)
  Y_pred <- as.numeric(model.response(mf_pred))
  rm(mf_orig)
  rm(mf_pred)
  col_name_keep <- colnames(X_orig)


  # Standardize to make estimation easier
  X_orig_use <- X_pred_use <- matrix(1, nrow = nrow(X_orig), ncol = ncol(X_orig))
  scale_cov  <- which(apply(X_pred, 2, function(x) sd(x, na.rm = TRUE)) > 0)

  if(all(X_orig[,1] == 1) == TRUE){
    with_intercept <- TRUE
    mean_X <- rep(0, ncol(X_pred))
    sd_X   <- rep(1, ncol(X_pred))
    for(j in scale_cov){
      X_pred_use[, j] <- scale(X_pred[, j])
      X_orig_use[, j] <- (X_orig[, j] - mean(X_pred[, j], na.rm = TRUE))/sd(X_pred[, j], na.rm = TRUE)
      mean_X[j] <- mean(X_pred[, j], na.rm = TRUE)
      sd_X[j] <- sd(X_pred[, j], na.rm = TRUE)
    }
    mean_X <- mean_X[-1] # remove the first element
    sd_X   <- sd_X[-1] # remove the first element
  }else{
    with_intercept <- FALSE

    sd_X <- rep(1, ncol(X_pred))
    for(j in scale_cov){
      X_pred_use[, j] <- X_pred[, j]/sd(X_pred[, j], na.rm = TRUE)
      X_orig_use[, j] <- X_orig[, j]/sd(X_pred[, j], na.rm = TRUE)
      sd_X[j] <- sd(X_pred[, j], na.rm = TRUE)
    }
  }
  colnames(X_orig_use) <- colnames(X_orig)
  colnames(X_pred_use) <- colnames(X_pred)

  # OLD
  # X_orig_use <- X_pred_use <- matrix(1, nrow = nrow(X_orig), ncol = ncol(X_orig))
  # # scale_cov  <- which(apply(X_pred, 2, function(x) sd(x, na.rm = TRUE)) > 0)
  # scale_cov  <- seq(from = 2, to = ncol(X_orig))
  # mean_X <- sd_X <- c()
  # for(j in scale_cov){
  #   X_pred_use[, j] <- scale(X_pred[, j])
  #   X_orig_use[, j] <- (X_orig[, j] - mean(X_pred[, j], na.rm = TRUE))/sd(X_pred[, j], na.rm = TRUE)
  #   mean_X[j-1] <- mean(X_pred[, j], na.rm = TRUE)
  #   sd_X[j-1] <- sd(X_pred[, j], na.rm = TRUE)
  # }
  # colnames(X_orig_use) <- colnames(X_orig)
  # colnames(X_pred_use) <- colnames(X_pred)

  rm(X_orig); rm(X_pred)

  # ##################################
  # Point Estimation
  # ##################################
  if(model == "felm"){
    # For Fixed Effects, we remove intercept
    X_orig_use <- X_orig_use[, -1, drop = FALSE]
    X_pred_use <- X_pred_use[, -1, drop = FALSE]

    # Fixed Effects
    adj_Y0 <- as.numeric(labeled_ind/sample_prob_use)*as.numeric(Y_pred - Y_orig)
    adj_Y0[labeled_ind == 0] <- 0
    adj_Y <- as.numeric(Y_pred) - adj_Y0
    adj_X0 <- (X_pred_use - X_orig_use) * as.numeric(labeled_ind/sample_prob_use)
    adj_X0[labeled_ind == 0, ] <- 0
    adj_X <- X_pred_use - adj_X0

    if(fixed_effect != "twoways"){
      fixed_effect_use <- demean_dsl(data_base = data_orig, adj_Y = adj_Y, adj_X = adj_X, index_use = index[1]) ## This part is correct.

      fe_Y <- fixed_effect_use$adj_Y_avg_exp
      fe_X <- fixed_effect_use$adj_X_avg_exp
    }else if(fixed_effect == "twoways"){
      fixed_effect_use_1 <- demean_dsl(data_base = data_orig, adj_Y = adj_Y, adj_X = adj_X, index_use = index[1])
      fixed_effect_use_2 <- demean_dsl(data_base = data_orig, adj_Y = adj_Y, adj_X = adj_X, index_use = index[2])

      fe_Y <- fixed_effect_use_1$adj_Y_avg_exp + fixed_effect_use_2$adj_Y_avg_exp - mean(adj_Y, na.rm = TRUE)
      fe_X <- t(t(fixed_effect_use_1$adj_X_avg_exp + fixed_effect_use_2$adj_X_avg_exp) - apply(adj_X, 2, mean, na.rm = TRUE))
    }
  }else{
    fe_Y <- fe_X <- NULL
  }

  inioptim <- rep(0, ncol(X_orig_use))

  est0 <- optim(par = inioptim,
                fn = dsl_general_moment,
                labeled_ind = labeled_ind,
                sample_prob_use = sample_prob_use,
                Y_orig = Y_orig,
                X_orig = X_orig_use,
                Y_pred = Y_pred,
                X_pred = X_pred_use,
                fe_Y = fe_Y,
                fe_X = fe_X,
                lambda = lambda,
                model = model,
                method = optim_method,
                hessian = FALSE,
                control = list(maxit = 5000),
                tuning_para = tuning_para)$par
  names(est0) <- colnames(X_orig_use)

  # ###################################################
  # (felm) we compute estimates for other paramters
  # ###################################################
  if(model == "felm"){
    if(fixed_effect == "oneway"){
      index_use <- index[1]
      adj_data <- fixed_effect_use$adj_data
      est_fe0   <- as.numeric(adj_data[, "adj_Y_avg"]) - as.matrix(adj_data[, colnames(X_orig_use), drop = FALSE])%*% est0
      formula_use <- paste0("~factor(",index_use, ")")
      X_fe <- model.matrixBayes(as.formula(formula_use), data = data_orig, drop.baseline = FALSE)
      colnames(X_fe) <- gsub(paste0("factor\\(",index_use, "\\)"), "", colnames(X_fe))
      est_fe <- est_fe0[match(colnames(X_fe), rownames(adj_data)), ]

      X_orig_use_exp <- cbind(X_orig_use, X_fe)
      X_pred_use_exp <- cbind(X_pred_use, X_fe)
      est0_exp <- c(est0, est_fe)
      rm(X_fe)
    }else if(fixed_effect == "twoways"){
      adj_data_1 <- fixed_effect_use_1$adj_data
      adj_data_2 <- fixed_effect_use_2$adj_data

      est_fe_base <- mean(adj_Y, na.rm = TRUE) - sum(est0* apply(adj_X, 2, mean, na.rm = TRUE))

      # alpha
      est_fe1_base   <- as.numeric(adj_data_1[, "adj_Y_avg"]) - as.matrix(adj_data_1[, colnames(X_orig_use), drop = FALSE])%*% est0 - est_fe_base
      formula_use1 <- paste0("~factor(",index[1], ")")
      X_fe1 <- model.matrixBayes(as.formula(formula_use1), data = data_orig, drop.baseline = FALSE)
      colnames(X_fe1) <- gsub(paste0("factor\\(",index[1], "\\)"), "", colnames(X_fe1))
      est_fe1_base <- est_fe1_base[match(colnames(X_fe1), rownames(adj_data_1)),]

      # gamma
      est_fe2_base   <- as.numeric(adj_data_2[, "adj_Y_avg"]) - as.matrix(adj_data_2[, colnames(X_orig_use), drop = FALSE])%*% est0
      formula_use2 <- paste0("~factor(",index[2], ")")
      X_fe2 <- model.matrixBayes(as.formula(formula_use2), data = data_orig, drop.baseline = FALSE)
      colnames(X_fe2) <- gsub(paste0("factor\\(",index[2], "\\)"), "", colnames(X_fe2))
      est_fe2_base <- est_fe2_base[match(colnames(X_fe2), rownames(adj_data_2)),]

      # we set gamma_0 = 0
      alpha_mean <- est_fe2_base[1]
      est_fe_1   <- est_fe1_base + alpha_mean
      est_fe_2   <- est_fe2_base - alpha_mean

      X_orig_use_exp <- cbind(X_orig_use, X_fe1, X_fe2[,-1, drop = FALSE])
      X_pred_use_exp <- cbind(X_pred_use, X_fe1, X_fe2[,-1, drop = FALSE])
      est0_exp <- c(est0, est_fe_1, est_fe_2[-1])
      rm(X_fe1); rm(X_fe2)
    }
  }else{
    est0_exp <- est0
    X_orig_use_exp <- X_orig_use
    X_pred_use_exp <- X_pred_use
  }
  rm(X_orig_use); rm(X_pred_use)

  # ################################
  # Variance Estimation
  # ################################
  n <- nrow(X_orig_use_exp)

  # Meat
  Meat_decomp <- dsl_general_moment_base_decomp(par = est0_exp,
                                                labeled_ind = labeled_ind,
                                                sample_prob_use = sample_prob_use,
                                                Y_orig = Y_orig,
                                                X_orig = X_orig_use_exp,
                                                Y_pred = Y_pred,
                                                X_pred = X_pred_use_exp,
                                                model  = model,
                                                clustered = clustered,
                                                cluster = data_orig$cluster___,
                                                tuning = tuning,
                                                tuning_para = tuning_para)
  # Meat   <- dsl_general_Meat(par = est0_exp, labeled_ind, sample_prob_use, Y_orig, X_orig_use_exp, Y_pred, X_pred_use_exp, model,
  #                            clustered, data_orig$cluster___)

  # Jacobian
  ## I do not incorporate the tuning parameter here because it does not affect asymptotically (2024/10/05)
  J <- dsl_general_Jacobian(par = est0_exp, labeled_ind, sample_prob_use, Y_orig, X_orig_use_exp, Y_pred, X_pred_use_exp, model)

  # Variance
  s_J <- solve(J)
  s_J_use <- s_J[1:length(est0), , drop = FALSE]
  rm(s_J)

  if(tuning == TRUE){

    # I think it is correct to pick tuning_para after standardization as we do here.(2024/10/05)

    tr_1 <- sum(diag(s_J_use %*% Meat_decomp$main_t1 %*% t(s_J_use)))
    tr_2 <- sum(diag(s_J_use %*% Meat_decomp$main_t2 %*% t(s_J_use)))

    tuning_para <- tr_2/(2*tr_1)

    tuning_para <- max(0, tuning_para)
    tuning_para <- min(1, tuning_para)
    return(tuning_para)

  }else if(tuning == FALSE){
    # Meat
    Meat <- Meat_decomp$main_1 + Meat_decomp$main_23

    V0 <- (s_J_use %*% Meat %*% t(s_J_use))/n

    # Variance of the target paramter
    if(model == "felm"){ # with intercept
      # Rescale-back the coefficients and V
      D_2 <- cbind(diag(1/sd_X, ncol = length(sd_X), nrow = length(sd_X))) # we don't need intercept
      D   <- rbind(D_2)
    }else{
      # Rescale-back the coefficients and V
      if(with_intercept == TRUE){
        D_1 <- c(1, - mean_X/sd_X)
        D_2 <- cbind(0, diag(1/sd_X, ncol = length(sd_X), nrow = length(sd_X)))
        D   <- rbind(D_1, D_2)
      }else{
        D_2 <- cbind(diag(1/sd_X, ncol = length(sd_X), nrow = length(sd_X)))
        D   <- rbind(D_2)
      }
    }

    # Estimates and Variance
    est <- as.numeric(D %*% est0)
    V   <- D %*% V0 %*% t(D)

    if(model == "felm"){
      col_name_keep <- col_name_keep[-1]
    }

    names(est) <- col_name_keep
    se <- sqrt(diag(V))
    names(se) <- col_name_keep

    out <- list("coefficients" = est, "standard_errors" = se, "vcov" = V,
                "Meat" = Meat, "Meat_decomp" = Meat_decomp, "J" = J, "D" = D, "vcov0" = V0)

  }


  return(out)
}
# ##############
# General
# ##############
dsl_general_moment <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_Y, fe_X, model, lambda, tuning_para){

  if(model == "lm"){
    m_dr   <- lm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)
  }else if(model == "logit"){
    m_dr   <- logit_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)
  }else if(model == "felm"){
    m_dr   <- felm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_Y, fe_X, tuning_para)
  }

  g <- apply(m_dr, 2, mean) # m_n(theta) in equation (6)
  g_out <- sum(g^2) + lambda*mean(par^2) # penalty

  return(g_out)
}

# dsl_general_Meat <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, clustered, cluster){
#
#   n <- nrow(X_orig)
#
#   if(model == "lm"){
#     m_dr   <- lm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
#   }else if(model == "logit"){
#     m_dr   <- logit_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
#   }else if(model == "felm"){ # we can use the same function as lm.
#     m_dr   <- lm_dsl_moment_base(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
#   }
#
#   if(clustered == FALSE){
#     Meat   <- (t(m_dr)%*% m_dr)/n
#   }else if(clustered == TRUE){
#     uniq_cluster <- sort(unique(cluster))
#     Meat_0 <- matrix(0, nrow = ncol(m_dr), ncol = ncol(m_dr))
#     for(z in 1:length(uniq_cluster)){
#       m_dr_cl <- m_dr[cluster == uniq_cluster[z], ,drop = FALSE]
#       matrix_ind <- matrix(1, nrow = nrow(m_dr_cl), ncol = nrow(m_dr_cl))
#       Meat_0 <- Meat_0 + t(m_dr_cl) %*% matrix_ind %*% m_dr_cl
#     }
#     Meat <- Meat_0/n
#   }
#
#   return(Meat)
# }

dsl_general_moment_base_decomp <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model, clustered, cluster, tuning, tuning_para){

  if(model == "lm"){
    m_orig <- lm_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)
    m_pred <- lm_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)

  }else if(model == "logit"){
    m_orig <- logit_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)
    m_pred <- logit_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)

  }else if(model == "felm"){ # we can use the same function as "lm"
    m_orig <- lm_dsl_moment_orig(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)
    m_pred <- lm_dsl_moment_pred(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, tuning_para)

    m_orig <- Matrix(m_orig, sparse = TRUE)
    m_pred <- Matrix(m_pred, sparse = TRUE)
  }

  n <- nrow(X_orig)

  # Clustering
  if(clustered == FALSE){
    # diag_1 <- Matrix(diag(labeled_ind/(sample_prob_use^2), ncol = length(labeled_ind), nrow = length(labeled_ind)), sparse = TRUE)
    # diag_2 <- Matrix(diag(labeled_ind/sample_prob_use, ncol = length(labeled_ind), nrow = length(labeled_ind)), sparse = TRUE)

    if(tuning == FALSE){
      diag_1 <- Diagonal(x = labeled_ind/(sample_prob_use^2))
      diag_2 <- Diagonal(x = labeled_ind/sample_prob_use)

      main_1 <- (t(m_orig - m_pred) %*% diag_1 %*% (m_orig - m_pred))/n
      main_2 <- (t(m_pred)%*% m_pred)/n
      main_30 <- t(m_pred)%*% diag_2 %*% (m_orig - m_pred)
      main_3 <- (main_30 + t(main_30))/n

      out   <- list("main_1" = main_1, "main_23" = main_2 + main_3)

    }else if(tuning == TRUE){
      diag_t1 <- Diagonal(x = c(1/sample_prob_use - 1))
      diag_t2 <- Diagonal(x = c(labeled_ind/(sample_prob_use^2) - labeled_ind/sample_prob_use))
      main_t1 <- (t(m_pred) %*% diag_t1 %*% (m_pred))/n  # E(1/pi m_pred*m_pred) - E(m_pred * m_pred)
      main_t2_0 <- (t(m_pred) %*% diag_t2 %*% (m_orig))/n  # E(1/pi m_pred*m_orig) - E(m_pred * m_orig)
      main_t2 <- main_t2_0 + t(main_t2_0)

      out <- list("main_t1" = main_t1, "main_t2" = main_t2)
    }

  }else if(clustered == TRUE){

    uniq_cluster <- sort(unique(cluster))

    if(tuning == FALSE){
      m_orig_m_pred <- (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
      m_pred_sum <- matrix(NA, ncol = ncol(m_pred), nrow = length(uniq_cluster))
      m_orig_m_pred_sum <- matrix(NA, ncol = ncol(m_orig_m_pred), nrow = length(uniq_cluster))
      for(z in 1:length(uniq_cluster)){
        m_pred_sum[z,1:ncol(m_pred)] <- colSums(m_pred[cluster == uniq_cluster[z], ,drop = FALSE])
        m_orig_m_pred_sum[z, 1:ncol(m_orig_m_pred)] <- colSums(m_orig_m_pred[cluster == uniq_cluster[z], ,drop = FALSE])
      }
      m_pred_sum <- Matrix(m_pred_sum, sparse = TRUE)
      m_orig_m_pred_sum <- Matrix(m_orig_m_pred_sum, sparse = TRUE)
      main_1  <- (t(m_orig_m_pred_sum)%*%m_orig_m_pred_sum)/n
      main_2  <- (t(m_pred_sum)%*%m_pred_sum)/n
      main_3_c0 <- (t(m_pred_sum)%*%m_orig_m_pred_sum)/n
      main_3  <- main_3_c0 + t(main_3_c0)

      out   <- list("main_1" = main_1, "main_23" = main_2 + main_3)

    }else if(tuning == TRUE){
      m_pred_t <- m_pred * as.numeric(sqrt(1/sample_prob_use - 1))
      m_pred_sum_t <- matrix(NA, ncol = ncol(m_pred), nrow = length(uniq_cluster))
      m_pred_t2 <- m_pred * as.numeric(c(labeled_ind/(sample_prob_use^2) - labeled_ind/sample_prob_use))
      m_pred_sum_t2 <- matrix(NA, ncol = ncol(m_pred), nrow = length(uniq_cluster))
      m_orig_sum_t <- matrix(NA, ncol = ncol(m_orig), nrow = length(uniq_cluster))
      for(z in 1:length(uniq_cluster)){
        m_pred_sum_t[z,1:ncol(m_pred_t)] <- colSums(m_pred_t[cluster == uniq_cluster[z], ,drop = FALSE])
        m_pred_sum_t2[z,1:ncol(m_pred_t)] <- colSums(m_pred_t2[cluster == uniq_cluster[z], ,drop = FALSE])
        m_orig_sum_t[z, 1:ncol(m_orig)] <- colSums(m_orig[cluster == uniq_cluster[z], ,drop = FALSE])
      }
      m_pred_sum_t <- Matrix(m_pred_sum_t, sparse = TRUE)
      main_t1  <- (t(m_pred_sum_t)%*%m_pred_sum_t)/n

      m_pred_sum_t2 <- Matrix(m_pred_sum_t2, sparse = TRUE)
      m_orig_sum_t  <- Matrix(m_orig_sum_t, sparse = TRUE)
      main_t2_0  <- (t(m_pred_sum_t2)%*%m_orig_sum_t)/n
      main_t2 <- main_t2_0 + t(main_t2_0)

      out <- list("main_t1" = main_t1, "main_t2" = main_t2)
    }

    # }else{
    # previous loop version
    #   uniq_cluster <- sort(unique(cluster))
    #   main_1_b <- main_2_b <- main_3_b <- matrix(0, nrow = ncol(m_pred), ncol = ncol(m_pred))
    #   for(z in 1:length(uniq_cluster)){
    #     m_orig_cl <- m_orig[cluster == uniq_cluster[z], ,drop = FALSE]
    #     m_pred_cl <- m_pred[cluster == uniq_cluster[z], ,drop = FALSE]
    #
    #     matrix_ind <- matrix(1, nrow = nrow(m_orig_cl), ncol = nrow(m_orig_cl))
    #
    #     m_orig_m_pred_cl <- (m_orig_cl - m_pred_cl) * as.numeric(labeled_ind/sample_prob_use)[cluster == uniq_cluster[z]]
    #
    #     if(model == "felm"){
    #       m_orig_m_pred_cl <- Matrix(m_orig_m_pred_cl, sparse = TRUE)
    #       matrix_ind <- Matrix(matrix_ind, sparse = TRUE)
    #     }
    #
    #     main_1_b <- main_1_b + (t(m_orig_m_pred_cl) %*% matrix_ind %*% (m_orig_m_pred_cl))
    #     main_2_b <- main_2_b + (t(m_pred_cl)%*% matrix_ind %*% m_pred_cl)
    #     main_3_b <- main_3_b + t(m_pred_cl)%*% matrix_ind %*% m_orig_m_pred_cl + t(m_orig_m_pred_cl)%*% matrix_ind %*% m_pred_cl
    #   }
    #   main_1 <- main_1_b/n
    #   main_2 <- main_2_b/n
    #   main_3 <- main_3_b/n
    # }
  }
  return(out)
}

dsl_general_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model){

  if(model == "lm"){
    J   <- lm_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model)
  }else if(model == "logit"){
    J   <- logit_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred)
  }else if(model == "felm"){
    # we can use the same function as "lm"
    J   <- lm_dsl_Jacobian(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model)
  }
  return(J)
}
