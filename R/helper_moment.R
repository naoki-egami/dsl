# ######
# lm
# ######
lm_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  m_orig <- X_orig * as.numeric(Y_orig - X_orig %*% par)

  m_orig[labeled_ind == 0, ] <- 0

  m_pred <- X_pred * as.numeric(Y_pred - X_pred %*% par)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

lm_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  m_orig <- X_orig * as.numeric(Y_orig - X_orig %*% par)
  m_orig[labeled_ind == 0, ] <- 0
  return(m_orig)
}

lm_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  m_pred <- X_pred * as.numeric(Y_pred - X_pred %*% par)
  return(m_pred)
}

lm_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, model){
  X_orig[labeled_ind == 0, ] <- 0

  if(model == "felm"){
    X_orig <- Matrix(X_orig, sparse = TRUE)
    X_pred <- Matrix(X_pred, sparse = TRUE)
  }
  diag_1 <- Diagonal(x = labeled_ind/sample_prob_use)
  diag_2 <- Diagonal(x = 1 - labeled_ind/sample_prob_use)

  # diag_1 <- diag(labeled_ind/sample_prob_use, ncol = length(labeled_ind), nrow = length(labeled_ind))
  # diag_2 <- diag(1 - labeled_ind/sample_prob_use, ncol = length(labeled_ind), nrow = length(labeled_ind))

  J <- (t(X_orig) %*% diag_1 %*% X_orig + t(X_pred) %*% diag_2 %*% X_pred)/nrow(X_orig)
  return(J)
}

# #####
# felm
# #####
felm_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred, fe_Y, fe_X){
  ## We just need to compute the alpha and gamma first.
  fe_use <- as.numeric(fe_Y) - as.matrix(fe_X) %*% as.numeric(par)

  m_orig <- X_orig * as.numeric(Y_orig - fe_use - X_orig %*% par)

  m_orig[labeled_ind == 0, ] <- 0

  m_pred <- X_pred * as.numeric(Y_pred - fe_use - X_pred %*% par)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

demean_dsl <- function(data_base, adj_Y, adj_X, index_use){
  data_base$uniq_id___ <- seq(1:nrow(data_base))
  adj_Y_avg  <- tapply(adj_Y, data_base[, index_use], mean, na.rm = TRUE)
  adj_X_avg0 <- lapply(1:ncol(adj_X), FUN = function(x) tapply(adj_X[, x], data_base[, index_use], mean, na.rm = TRUE))
  adj_X_avg  <- do.call("cbind", adj_X_avg0)
  colnames(adj_X_avg) <- colnames(adj_X)
  adj_data   <- data.frame(cbind(adj_Y_avg, adj_X_avg))
  adj_data$index___ <- rownames(adj_data)
  fixed_effect_use <- merge(data_base[, c(index_use, "uniq_id___"), drop = FALSE], adj_data, by.x = index_use, by.y = "index___", all.x = TRUE, sort = FALSE)
  fixed_effect_use <- fixed_effect_use[order(fixed_effect_use$uniq_id___, decreasing = FALSE), , drop = FALSE]
  rm(data_base)
  out <- list("adj_Y_avg_exp" = fixed_effect_use[, "adj_Y_avg"], "adj_X_avg_exp" = fixed_effect_use[, colnames(adj_X), drop = FALSE],
              "adj_data" = adj_data)
  return(out)
}

# #####
# logit
# #####
logit_dsl_moment_base <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  inv_logit <- 1/(1+exp(-X_orig %*% par))
  m_orig <- X_orig * as.numeric(Y_orig - inv_logit)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y

  inv_logit_pred <- 1/(1+exp(-X_pred %*% par))
  m_pred <- X_pred * as.numeric(Y_pred - inv_logit_pred)
  m_dr   <- m_pred + (m_orig - m_pred) * as.numeric(labeled_ind/sample_prob_use)
  return(m_dr)
}

logit_dsl_moment_orig <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  inv_logit <- 1/(1+exp(-X_orig %*% par))
  m_orig <- X_orig * as.numeric(Y_orig - inv_logit)
  m_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  return(m_orig)
}

logit_dsl_moment_pred <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){
  inv_logit_pred <- 1/(1+exp(-X_pred %*% par))
  m_pred <- X_pred * as.numeric(Y_pred - inv_logit_pred)
  return(m_pred)
}

logit_dsl_Jacobian <- function(par, labeled_ind, sample_prob_use, Y_orig, X_orig, Y_pred, X_pred){

  #
  # inv_logit_pred <- 1/(1+exp(-X_pred%*%par))
  # inv_logit_pred2 <- inv_logit_pred/(1 + inv_logit_pred)^2
  inv_logit_pred <- exp(-X_pred%*%par)
  inv_logit_pred2 <- inv_logit_pred/(1 + inv_logit_pred)^2
  diag_pred2 <- Diagonal(x = inv_logit_pred2)
  diag_pred2_R <- Diagonal(x = inv_logit_pred2*labeled_ind/sample_prob_use)

  grad_pred <- (t(X_pred) %*% diag_pred2 %*% X_pred)/nrow(X_pred)
  grad_pred_R <- (t(X_pred) %*% diag_pred2_R %*% X_pred)/nrow(X_pred)

  # original
  # inv_logit_orig  <- 1/(1+exp(-X_orig%*%par))
  # # inv_logit_orig2 <- inv_logit_pred/(1 + inv_logit_pred)^2
  # inv_logit_orig2 <- inv_logit_orig/(1 + inv_logit_orig)^2
  inv_logit_orig  <- exp(-X_orig%*%par)
  inv_logit_orig2 <- inv_logit_orig/(1 + inv_logit_orig)^2
  inv_logit_orig2[labeled_ind == 0] <- 0  # r/pi * Y
  X_orig[labeled_ind == 0, ] <- 0  # r/pi * Y
  diag_orig2_R <- Diagonal(x = inv_logit_orig2 *labeled_ind/sample_prob_use)
  grad_orig <- (t(X_orig) %*% diag_orig2_R %*% X_orig)/nrow(X_orig)

  grad_main <- grad_pred + grad_orig - grad_pred_R

  out <- grad_main
  return(out)
}
