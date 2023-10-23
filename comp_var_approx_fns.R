# functions for computing variational approximation
# performing update steps for fitting variational and model parameters

# version of parLapply that will revert to regular lapply if cl is NULL
lapply_parLapply <- function(cl, X, FUN){
  if(is.null(cl)){
    lapply(X, FUN)
  }else{
    parLapply(cl, X, FUN)
  }
}

# version of parLapplyLB that will revert to regular lapply if cl is NULL
lapply_parLapplyLB <- function(cl, X, FUN){
  if(is.null(cl)){
    lapply(X, FUN)
  }else{
    parLapplyLB(cl, X, FUN)
  }
}

### HELPER FUNCTION
scaled_sum <- function(rate, vec=NULL, len=NULL){
  if(is.null(len)){
    len <- length(vec)
  }

  if (is.null(vec)){
    if(rate == 1){
      return(len)
    } else {
      return((1 - rate ^ len) / (1 - rate))
    }
  } else {
    if(rate == 1){
      return(sum(vec))
    } else {
      return(sum(vec * rate ^ ((len - 1):0)))
    }
  }
}

# apply cumulative scaling more efficiently to cols of matrix
scaled_sum_mat <- function(rate, mat){
  if(rate != 1 & nrow(mat) > 1){
    for (row_k in 2:nrow(mat)){
      mat[row_k, ] <- rate * mat[row_k - 1, ] + mat[row_k, ]
    }
  }
  return(mat)
}
# return scaled_sum_mat multiplied by vector 
scaled_sum_mat_mult <- function(rate, mat, vec){
  n <- nrow(mat)
  out_vec <- numeric(n)
  for (k in 1:n){
    if (k > 1 & rate != 1){
      mat[k, ] <- rate * mat[k - 1, ] + mat[k, ]
    }
    out_vec[k] <- sum(vec * mat[k, ])
  }
  return(out_vec)
}

# get full alpha vector and gradients w.r.t. nu1, nu2
# T_: inter spike interval length, "spike time" following reset
get_alphas <- function(T_, nu1, nu2, delta, retGrad = TRUE){
  
  no_nu1 <- exp(-nu2 * ((T_ - 1):0)) - 
    ((delta + 1) * scaled_sum(rate = delta * exp(-nu2), len = T_) - 1) /
      ((delta + 1) * scaled_sum(rate = delta, len = T_) - 1)
  
  alphasOut <- nu1 * no_nu1
  
  gradOut <- NULL
  
  if(retGrad){
    gradOut <- matrix(nrow = T_, ncol = 2)
    colnames(gradOut) <- c("nu1", "nu2")
    
    gradOut[, 1] <- no_nu1
    gradOut[, 2] <- -nu1 * ((T_ - 1):0) * exp(-nu2 * ((T_ - 1):0)) +
      nu1 * (delta + 1) * scaled_sum(rate = exp(-nu2) * delta, vec = (T_ - 1):0) /
        ((delta + 1) * scaled_sum(rate = delta, len = T_) - 1)
  }
  
  return(
    list(
      alphas = alphasOut,
      alphas_grad_mat = gradOut
    )
  )
}

# reparameterization of phi so we can use unconstrained gradient steps
get_reparam_phi <- function(phi, ratio_bound_for_reparam = .25){
  ratio_bound_for_reparam = ifelse(is.null(ratio_bound_for_reparam), .25, ratio_bound_for_reparam)
  
  reparam_phi_names <- c("log_vartheta", "log_nu1", "log_nu2", "logit_u1", "log_u2", "logit_beta1", "log_beta2")
  phi_reparam <- numeric(7)
  names(phi_reparam) <- reparam_phi_names
  phi_reparam[c("log_vartheta", "log_nu1", "log_nu2", "log_beta2")] <-
    log(phi[c("vartheta", "nu1", "nu2", "beta2")])
  phi_reparam["logit_beta1"] <- gtools::logit(phi["beta1"])
  u_1 <- (phi["sigma"] ^ 2 / phi["tau"] ^ 2) * ratio_bound_for_reparam
  phi_reparam["logit_u1"] <- gtools::logit(u_1)
  phi_reparam["log_u2"] <- log(phi["tau"])
  return(phi_reparam)
}

# inverse of reparameterization, to recover phi
get_phi_from_reparam <- function(phi_reparam, ratio_bound_for_reparam = .25){
  ratio_bound_for_reparam = ifelse(is.null(ratio_bound_for_reparam), .25, ratio_bound_for_reparam)
  
  phi_names <- c("vartheta", "nu1", "nu2", "sigma", "tau", "beta1", "beta2")
  phi <- numeric(7)
  names(phi) <- phi_names
  phi[c("vartheta", "nu1", "nu2", "beta2")] <-
    exp(phi_reparam[c("log_vartheta", "log_nu1", "log_nu2", "log_beta2")])
  phi["beta1"] <- gtools::inv.logit(phi_reparam["logit_beta1"])
  
  u_1 <- gtools::inv.logit(phi_reparam["logit_u1"])
  u_2 <- exp(phi_reparam[c("log_u2")])
  
  phi["sigma"] <- u_2 * sqrt(u_1 / ratio_bound_for_reparam)
  phi["tau"] <- u_2
  
  return(phi)
}

# get derivatives wrt reparameterized variables
get_grad_reparam_phi <- function(grad_phi, phi, ratio_bound_for_reparam = .25){
  ratio_bound_for_reparam = ifelse(is.null(ratio_bound_for_reparam), .25, ratio_bound_for_reparam)
  
  reparam_phi_names <- c("log_vartheta", "log_nu1", "log_nu2", "logit_u1", "log_u2", "logit_beta1", "log_beta2")
  grad_reparam <- numeric(length(reparam_phi_names))
  names(grad_reparam) <- reparam_phi_names
  
  grad_reparam[c("log_vartheta", "log_nu1", "log_nu2", "log_beta2")] <-
    grad_phi[c("vartheta", "nu1", "nu2", "beta2")] * phi[c("vartheta", "nu1", "nu2", "beta2")]
  
  grad_reparam["logit_beta1"] <- grad_phi["beta1"] * phi["beta1"] * (1 - phi["beta1"])
  
  grad_reparam["logit_u1"] <- grad_phi["sigma"] * (phi["sigma"] / 2) *
    (1 -  ratio_bound_for_reparam * phi["sigma"] ^ 2 / phi["tau"] ^ 2)
  
  grad_reparam["log_u2"] <- grad_phi["tau"] * phi["tau"] +
    grad_phi["sigma"] * phi["sigma"]
  
  return(grad_reparam)
}

# functions for moving between variables and reparamerized variables
get_reparam_from_variable <- function(var_in, var_name, param_list){
  if(var_name == "phi"){
    if(is.null(dim(var_in))){
      get_reparam_phi(var_in, ratio_bound_for_reparam=param_list$ratio_bound_for_reparam)
    }else{
      sapply(
        1:ncol(var_in),
        function(node_i){
          get_reparam_phi(
            var_in[, node_i], ratio_bound_for_reparam=param_list$ratio_bound_for_reparam
          )
        }
      )
    }
  }else if(var_name == "W"){
    if (param_list$restrict_W){
      gtools::logit(var_in, min=param_list$W_lb, max=param_list$W_ub)
    } else {
      var_in
    }
  }else if(var_name == "eta"){
    gtools::logit(var_in, min=param_list$eta_lb, max=param_list$eta_ub)
  }
}

get_variable_from_reparam <- function(var_in, var_name, param_list){
  if(var_name == "phi"){
    if(is.null(dim(var_in))){
      get_phi_from_reparam(var_in, ratio_bound_for_reparam=param_list$ratio_bound_for_reparam)
    }else{
      sapply(
        1:ncol(var_in),
        function(node_i){
          get_phi_from_reparam(
            var_in[, node_i], ratio_bound_for_reparam=param_list$ratio_bound_for_reparam
          )
        }
      )
    }
  }else if(var_name == "W"){
    if (param_list$restrict_W){
      gtools::inv.logit(var_in, min=param_list$W_lb, max=param_list$W_ub)
    } else {
      var_in
    }
  }else if(var_name == "eta"){
    gtools::inv.logit(var_in, min=param_list$eta_lb, max=param_list$eta_ub)
  }
}

# get data for indexed inter-spike interval
# use interval index counter
# X is the full spike train data matrix, 0 or 1s
get_interval_data <- function(interval_index, spike_intervals, binary_spike_train){
  spike_interval_data <- spike_intervals %>% filter(index == interval_index)
  return(
    list(
      i = spike_interval_data$i,
      T_ = spike_interval_data$T_,
      X = binary_spike_train[spike_interval_data$t_start:(spike_interval_data$t_end - 1), -spike_interval_data$i, drop = FALSE]
    )
  )
}

#### helper functions for computing entries of the Cholesky decompition matrix R that defines the covariance of q
# compute diagonal entries of R, unmultiplied by sigma
get_r_t <- function(t, T_, phi, grads=TRUE){
  t_T_ratio <- (T_ - t) / (T_ - 1)
  exp_term_numerator <- exp(-phi["beta2"] * t_T_ratio)
  numerator <- 1 - phi["beta1"] * exp_term_numerator
  denominator <- (1 - phi["beta1"] * exp(-phi["beta2"]))
  
  r_t <- unname(numerator / denominator)
  
  if (grads) {
    grad_beta1 <- (exp(-phi["beta2"]) - exp_term_numerator) / denominator ^ 2
    grad_beta2 <- 
      phi["beta1"] * (t_T_ratio * exp_term_numerator - exp(-phi["beta2"]) * (1 - (1 - t_T_ratio) * phi["beta1"] * exp_term_numerator)) /
      denominator ^ 2
    
    return(cbind("r_t" = r_t, "beta1" = unname(grad_beta1), "beta2" = unname(grad_beta2)))
  } else {
    return(cbind("r_t" = r_t))
  }
}

# compute unscaled off-diagonal (below diagonal) entries of R, unmultiplied by sigma
getL <- function(T_, phi, c_vec, grads=TRUE){
  r_t_out <- get_r_t(t = 1:T_, T_ = T_, phi = phi, grads = grads)
  
  retList <- list()
  retList$L <- -((1 / c_vec) %*% t(c_vec)) %*% diag(r_t_out[, "r_t"] / c(T_ - 1:(T_ - 1), 1)) * lower.tri(diag(T_))
  if (grads){
    for (beta_i in c("beta1", "beta2")){
      retList[[beta_i]] <- -((1 / c_vec) %*% t(c_vec)) %*% diag(r_t_out[, beta_i] / c(T_ - 1:(T_ - 1), 1)) * lower.tri(diag(T_))
    }
  }
  return(retList)
}

# compute the scale for the off-diagonal elements of R, a function of
# beta1, beta2, tau, sigma, delta and T_
# note that we must have phi["tau"] / phi["sigma"] >= c_vec[T_] * r_T[1, "r_t"],
# this is ensured by the reparameterization of phi
get_scale <- function(T_, phi, c_vec, grads=TRUE, return_L_out = TRUE){
  retVals <- list()
  
  # first compute the scalar "a", from the quadratic form of L and c_vec
  L_out <- getL(T_ = T_, c_vec=c_vec, phi = phi, grads = grads)
  a <- quad.form(L_out$L, c_vec, chol = T)
  sqrt_a <- sqrt(a)
  
  # gradients of a wrt beta1 and beta2
  if(grads){
    grads_a <- list()
    for (beta_i in c("beta1", "beta2")){
      grads_a[beta_i] <-
        quad.form(L_out$L %*% t(L_out[[beta_i]]) + L_out[[beta_i]] %*% t(L_out$L), c_vec, chol = F) / (2 * sqrt_a)
    }
  }
  
  # now we compute the scale with the quadratic formula
  r_T <- get_r_t(t = T_, T_ = T_, phi = phi, grads = grads)
  sqrt_term <- sqrt((phi["tau"] / phi["sigma"]) ^ 2 - (c_vec[T_] * r_T[1, "r_t"]) ^ 2)
  retVals$scale <- 1 - sqrt_term / sqrt_a
  
  if(grads){
    for (beta_i in c("beta1", "beta2")){
      retVals[beta_i] <-
        c_vec[T_] ^ 2 * r_T[1, "r_t"] * r_T[1, beta_i] / (sqrt_a * sqrt_term) + 
        sqrt_term * grads_a[[beta_i]] / a
    }
    
    retVals$sigma <- (phi["tau"] ^ 2 / phi["sigma"] ^ 3) / (sqrt_a * sqrt_term)
    retVals$tau <- -1 * (phi["tau"] / phi["sigma"] ^ 2) / (sqrt_a * sqrt_term)
  }
  
  if (return_L_out){
    retVals$L_out <- L_out
  }
  
  return(retVals)
}

# compute the Cholesky decompition matrix R, using helper functions above
# return additional pieces including determinant
# and quadratic forms giving variances of V_t
# for getting R for simulation, call "get_R(T_ = T_, phi = phi, delta=delta, c_vec=c_vec, grads=FALSE, ret_names="R")"
get_R <- function(T_, phi, delta, c_vec, grads=TRUE, get_s_vec_out = T, ret_names = c("det_out", "variation_out")){
  retList <- list()
  
  # take care of trivial case where T_ = 1
  if (T_ == 1){
    retList$R <- matrix(phi["sigma"], nrow = 1, ncol = 1)
    retList$det_out <- list(halfLogDet = log(phi["sigma"]))
    if(grads){
      for (beta_i in c("beta1", "beta2")){
        retList[[beta_i]] <- matrix(0, nrow = 1, ncol = 1)
      }
      retList$sigma <- matrix(1, nrow = 1, ncol = 1)
      retList$tau <- matrix(0, nrow = 1, ncol = 1)
      retList$det_out$sigma <- unname(T_ / phi["sigma"])
      retList$det_out$tau <- 0
      for (beta_i in c("beta1", "beta2")){
        retList$det_out[[beta_i]] <- 0
      }
    }
    if (get_s_vec_out){
      retList$variation_out <-
        cbind("t" = 1, "s_vec_sqrd" = phi["sigma"] ^ 2, "s_vec" = phi["sigma"], "Sigma_diag" = phi["sigma"] ^ 2)
      if(grads){
        Sigma_diag_grad <- sapply(
          c("beta1", "beta2", "sigma", "tau"),
          function(curVar){
            Rg_t <- retList$R %*% t(retList[[curVar]])
            return(diag(Rg_t + t(Rg_t)))
          }
        ) %>% matrix(nrow = 1)
        colnames(Sigma_diag_grad) <- paste0("Sigma_diag_grad_", c("beta1", "beta2", "sigma", "tau"))
        retList$variation_out <- cbind(retList$variation_out, Sigma_diag_grad)
      }
    }
    return(retList[ret_names])
  }
  
  r_t_out <- get_r_t(t = 1:T_, T_ = T_, phi = phi, grads = grads)
  scale_out <- get_scale(T_ = T_, phi = phi, c_vec=c_vec, grads = grads, return_L_out = TRUE)
  
  # first fill in off diagonal starting from L_out
  retList$R <- phi["sigma"] * scale_out$scale * scale_out$L_out$L
  # then add the diagonal entries
  diag(retList$R) <- phi["sigma"] * r_t_out[, "r_t"]
  
  if (grads){
    for (beta_i in c("beta1", "beta2")){
      # first fill in off diagonal starting from L_out
      retList[[beta_i]] <- phi["sigma"] *
        (scale_out$scale * scale_out$L_out[[beta_i]] + scale_out[[beta_i]] * scale_out$L_out$L)
      # then add the diagonal entries
      diag(retList[[beta_i]]) <- phi["sigma"] * r_t_out[, beta_i]
    }
    
    ### next gradient mat wrt sigma
    # fill in off diagonal starting from L_out
    retList$sigma <- (scale_out$scale + scale_out[["sigma"]] * phi["sigma"]) * scale_out$L_out$L
    # then add the diagonal entries
    diag(retList$sigma) <- r_t_out[, "r_t"]
    
    ### finally gradient mat wrt tau (note zeros along diagonal)
    retList$tau <- scale_out[["tau"]] * phi["sigma"] * scale_out$L_out$L
  }
  
  # compute the 1/2 log determinant of Sigma = R %*% R^T
  retList$det_out <- list()
  retList$det_out$halfLogDet <- unname(T_ * log(phi["sigma"]) + sum(log(r_t_out[, "r_t"])))
  if(grads){
    retList$det_out$sigma <- unname(T_ / phi["sigma"])
    retList$det_out$tau <- 0
    for (beta_i in c("beta1", "beta2")){
      retList$det_out[[beta_i]] <- sum(r_t_out[, beta_i] / r_t_out[, "r_t"])
    }
  }
  
  # get sd vec of v_t, t=1 to t=T_ under q_theta
  if (get_s_vec_out){
    
    # initialize stored values and matrix for storing output
    if(grads){
      grad_vars <- c("beta1", "beta2", "sigma", "tau")
      # we will be storing multiplied matrices, updating recursively
      # taking advantage of the fact that R (along with its gradients) is lower-triangular
      
      # initialize the stored lists
      var_grads_cur <- numeric(4)
      names(var_grads_cur) <- grad_vars
      
      s_vec_out <- matrix(0, nrow = T_, ncol = 12)
      colnames(s_vec_out) <- c("t", "s_vec_sqrd", "s_vec", "Sigma_diag",
                               paste0("s_vec_grad_", grad_vars),
                               paste0("Sigma_diag_grad_", grad_vars))
    } else {
      s_vec_out <- matrix(0, nrow = T_, ncol = 4)
      colnames(s_vec_out) <- c("t", "s_vec_sqrd", "s_vec", "Sigma_diag")
    }
    
    for (t in 1:T_){
      s_vec_out[t, "t"] <- t
      
      # compute the Sigma_diag entries
      s_vec_out[t, "Sigma_diag"] <- sum(retList$R[t, 1:t] ^ 2) # t bits
      if(grads){
        for (curVar in grad_vars){
          s_vec_out[t, paste0("Sigma_diag_grad_", curVar)] <- 2 * sum(retList$R[t, 1:t] * retList[[curVar]][t, 1:t])
        }
      }
      
      # start with base case t = 1
      if (t == 1) {
        var <- retList$R[1, 1] ^ 2
        sd <- sqrt(var)
        s_vec_out[t, "s_vec_sqrd"] <- var
        s_vec_out[t, "s_vec"] <- sd
        
        if(grads){
          # initialize stored quadratic forms matrix multiplication values
          for (curVar in grad_vars){
            # calculate and store base case output value
            var_grads_cur[[curVar]] <- 2 * retList$R[1, 1] * retList[[curVar]][1, 1]
            s_vec_out[t, paste0("s_vec_grad_", curVar)] <- var_grads_cur[[curVar]] / (2 * sd)
          }
        }
        
      } else {
        var <- (delta ^ 2) * s_vec_out[t - 1, "s_vec_sqrd"] +
          2 * sum(retList$R[1:(t - 1), 1:(t - 1), drop = F] %*% t(retList$R[t, 1:(t - 1), drop = F]) * delta ^ ((t - 1):1)) + 
          s_vec_out[t, "Sigma_diag"]
        
        sd <- sqrt(var)
        s_vec_out[t, "s_vec_sqrd"] <- var
        s_vec_out[t, "s_vec"] <- sd
        
        if(grads){
          for (curVar in grad_vars){
            # update stored matrix multiplication
            row_v <- retList$R[t, 1:(t - 1), drop = F] %*% t(retList[[curVar]][1:(t - 1), 1:(t - 1), drop = F])
            col_v <- retList$R[1:(t - 1), 1:(t - 1), drop = F] %*% t(retList[[curVar]][t, 1:(t - 1), drop = F])
            
            var_grads_cur[[curVar]] <-
              (delta ^ 2) * var_grads_cur[[curVar]] +
              2 * sum((as.numeric(row_v) + as.numeric(col_v)) * delta ^ ((t - 1):1)) +
              s_vec_out[t, paste0("Sigma_diag_grad_", curVar)]
            
            s_vec_out[t, paste0("s_vec_grad_", curVar)] <- var_grads_cur[[curVar]] / (2 * sd)
          }
        }
      }
    }
    retList$variation_out <- s_vec_out
  }
  
  return(retList[ret_names])
}

### MAIN FUNCTION
# function that calculates objective for single inter spike interval
# aand gradients with respect to W, phi, eta
# inputs: 
#     W_in: matrix edge weights into spiking node (vector of length d = N - 1, where N size of network)
#     phi: named vector with entries "sigma", "vartheta", "tau_1", "nu1", "nu2"
#     eta: scalar constant input
#     T_: inter spike interval length
#     X: T_ by d matrix of input spike indicators
#     grads: specifies which set of variable to calculate gradients wrt (change to all or none...)
obj_grad_single_spike <- function(
  W_in, phi, eta,
  T_, X,
  delta, K=50, piecewise_approx_list = NULL,
  grads = TRUE,
  retFittedVals = FALSE,
  R_out = NULL  #option to provide R_out pre-calculated
){
  retList <- list()
  d <- length(W_in)
  
  # get cumulative/scaled sum of input *unmoderated* by W_in
  # in neural model (LIF), need scaling / discount function not just 'cumsum'
  # the t^th row of X_scaled is the vector of cumulative (and scaled) binary input spikes,
  # t steps into the (selected, observed) interspike interval
  if (grads){
    scaled_cum_spike_mat <- rbind(0, scaled_sum_mat(delta, X))
    get_R_1_mat_fn <- get_R_1_mat
    get_R_2_mat_fn <- get_R_2_mat
  }
  else {
    get_R_1_mat_fn <- get_R_1_mat_6
    get_R_2_mat_fn <- get_R_2_mat_6
  }
  
  # cumulative/scaled sum of input *moderated* by W_in
  scaled_cum_arrivals <- c(0, scaled_sum_mat_mult(delta, X, W_in))
  
  # sum_{k=1}^r delta^(r-k) vec for r=1 to T_
  scaled_sums_store <- (1 - delta ^ (1:T_)) / (1 - delta)
  
  # and sum_{k=1}^r delta^2(r-k) vec for r=1 to T_
  scaled_sums_store_2 <- (1 - delta ^ (2 * (1:T_))) / (1 - delta ^ 2)
  
  # total nonrandom input
  non_z_cum_input_vec <- scaled_cum_arrivals + eta * scaled_sums_store
  
  # get alphas vector
  # get gradients wrt nu1 and nu2 if grads == "phi"
  alphasOut <- get_alphas(T_, phi["nu1"], phi["nu2"], delta, retGrad = grads)
  
  # save key scaled sums
  s0 <- scaled_sum(rate = delta, len = T_)
  s1 <- scaled_sum(rate = delta, len = T_) ^ 2
  s2 <- scaled_sum(rate = delta ^ 2, len = T_)
  
  ### need to calculate approximate contribution from log p(y_t | v_t) all t = 1,...,T_
  ### working towards this we compute vectors of means and variances of z_t and v_t for t = 1,2,...,T_
  # expectations E z_t, t=1 to t=T_ under q_phi
  c_mu <- (2 * phi["vartheta"] - non_z_cum_input_vec[T_] - non_z_cum_input_vec[T_ - 1]) /
    ((1 + delta) * s2 - 1)
  mus <- c_mu * delta ^ (T_:1) + alphasOut$alphas
  # related quanities
  sum_mus <- sum(mus)
  sum_mus_sqrd <- sum(mus ^ 2)
  # next calculate mean vec E v_t, t=1 to t=T_ under q_phi
  # can check that vartheta == .5 * (m_vec[T_] + m_vec[T_ - 1])
  scaled_cumsums_mus <- vapply(1:T_, function(ell) scaled_sum(rate = delta, vec = mus[1:ell]), 0)
  m_vec <- non_z_cum_input_vec + scaled_cumsums_mus
  
  # Next we get needed calculated quantities from the Cholesky decompition matrix R that defines
  # the covariance matrix Sigma = R %*% t(R) under q_phi
  
  # set c_vec s.t. var(V* = .5 * (V_T + V_{T-1})) = tau ^ 2
  c_vec <- .5 * (delta ^ ((T_ - 1):0) + c(delta ^ ((T_ - 2):0), 0))
  
  # if just interested in fitted values, return those now (and exit function)
  if(retFittedVals){
    R_out = get_R(T_ = T_, phi = phi, delta = delta, c_vec = c_vec,
                  grads = FALSE, get_s_vec_out = TRUE, ret_names = c("R", "variation_out"))
    return(
      list(
        mus = mus,
        R = R_out$R,
        non_z_cum_input_vec = non_z_cum_input_vec,
        m_vec = m_vec,
        s_vec = R_out$variation_out[, "s_vec"],
        R_out = R_out
      )
    )
  }
  
  # otherwise get determinant and s_vec, with derivatives
  if (is.null(R_out)){
    R_out <- get_R(T_ = T_, phi = phi, delta = delta, c_vec = c_vec, grads = grads, get_s_vec_out = TRUE)
  }
  # sds of v_t for t = 1,2,...,T_
  s_vec <- R_out$variation_out[, "s_vec"]
  # variances of z_t for t = 1,2,...,T_
  z_vars <- R_out$variation_out[, "Sigma_diag"]
  
  ## store ingredients for calculating E_xk_vec_mat
  knots <- lapply(piecewise_approx_list, "[[", "interval") %>% unlist %>% unique
  R1_mat <- get_R_1_mat_fn(mu=m_vec, sigma=s_vec)
  R2_mat_list <- lapply(knots, function(knot) get_R_2_mat_fn(x=knot, mu=m_vec, sigma=s_vec))
  pn_list <- lapply(knots, function(knot) pnorm(knot, m_vec, s_vec))
  dn_list <- lapply(knots, function(knot) dnorm(knot, m_vec, s_vec))
  
  # NOW calculate sum E[h(V_t)] from t = 1 to t = T_ and its gradient wrt W_in
  # to do this, for each segment of the piecewise approximation
  # get segment's interval (a,b) and the polynomial coeffiecients of h on the segment
  # calculate the E[V_t^k * ind_[a,b](x)] for t=1,...,T_ and k = 0,...,7 and store in a 7 x T_ matrix 
  # get objective and gradient contributions for the given segment
  # add up over all the segements
  
  E_hV_w_grads_segments <- lapply(
    1:length(piecewise_approx_list),
    function(segment_index){
      
      segmentList <- list()
      segment <- piecewise_approx_list[[segment_index]]
      
      # much more efficient vectorized computation
      E_xk_vec_mat <- t(2 * (pn_list[[segment_index + 1]] - pn_list[[segment_index]]) * R1_mat -
                          dn_list[[segment_index + 1]] * R2_mat_list[[segment_index + 1]] +
                          dn_list[[segment_index]] * R2_mat_list[[segment_index]])
      
      # objecitve contribution (negative)
      segmentList$neg_obj_contrib <- sum(E_xk_vec_mat[1:length(segment$coefs), ] * segment$coefs)
      
      if (grads){
        EpX_vec <- colSums(E_xk_vec_mat[1:length(segment$coefs), ] * segment$coefs)
        EXpX_vec <- colSums(E_xk_vec_mat[2:(length(segment$coefs) + 1), ] * segment$coefs)
        EX2pX_vec <- colSums(E_xk_vec_mat[3:(length(segment$coefs) + 2), ] * segment$coefs)
        
        segmentList$integrals_m <- (EXpX_vec - m_vec * EpX_vec) / (s_vec ^ 2)
        segmentList$integrals_s <-
          (EX2pX_vec / s_vec ^ 2 - EXpX_vec * 2 * m_vec / s_vec ^ 2 + EpX_vec * (m_vec ^ 2 / s_vec ^ 2  - 1)) / s_vec
      }
      
      return(segmentList)
    }
  )
  
  # sum over segments of piecewise approximation h of f
  E_hV_contrib <- list()
  for (contrib_name in names(E_hV_w_grads_segments[[1]])){
    E_hV_contrib[[contrib_name]] <- Reduce("+", lapply(E_hV_w_grads_segments, "[[", contrib_name))
  }
  
  # calculate objective in single step
  obj <- 
      # contribution from E_{q(z)} sum_t log p(z_t)
    -T_ * (log(phi["sigma"]) + log(2 * pi) / 2) -
    (sum(z_vars) + sum_mus_sqrd) / (2 * phi["sigma"] ^ 2) -
      # subtract sum of contributions from E[h(V_t)]
    E_hV_contrib$neg_obj_contrib +
      # extra component from E[log p(y[T] = 1 | V_T)]: kappa * (E_q v_T - 1) 
    K * (2 * delta * phi["vartheta"] + eta + mus[T_] + sum(W_in * X[T_ - 1, ])) / (1 + delta) - K +
      # entropy of q = .5 * log det Sigma
    R_out$det_out$halfLogDet
  
  retList$obj <- unname(obj)
  
  # part of objective fuction that DOES NOT depend on input from network
  obj_component0 <- 
    -T_ * (log(phi["sigma"]) + log(2 * pi) / 2) -
    sum(z_vars) / (2 * phi["sigma"] ^ 2) +
    K * (2 * delta * phi["vartheta"] + eta) / (1 + delta) - K +
    R_out$det_out$halfLogDet
  
  # part of objective fuction that DOES depend on input from network
  obj_component1 <- 
    -sum_mus_sqrd / (2 * phi["sigma"] ^ 2) -
    E_hV_contrib$neg_obj_contrib + K * (mus[T_] + sum(W_in * X[T_ - 1, ])) / (1 + delta)
  
  retList$obj_component0 <- unname(obj_component0)
  retList$obj_component1 <- unname(obj_component1)
  
  # calculagte gradients to return (if any)
  if(grads){
    
    #### phi GRADIENT FIRST
    # initialize empty gradient vector
    grad_phi <- numeric(7)
    names(grad_phi) <- names(phi)
    
    ## add gradients contribution of E_{q(z)} sum_t log p(z_t)
    grad_phi["vartheta"] <- ( -1 / (phi["sigma"] ^ 2) ) *
      2 * delta * scaled_cumsums_mus[T_] / ((1 + delta) * s2 - 1)
    grad_phi[c("nu1", "nu2")] <- ( -1 / (phi["sigma"] ^ 2) ) *
      colSums(alphasOut$alphas_grad_mat * mus)
    
    for (curVar in c("tau", "beta1", "beta2")){
      grad_phi[curVar] <- -.5 * sum(R_out$variation_out[, paste0("Sigma_diag_grad_", curVar)]) / phi["sigma"] ^ 2
    }
    grad_phi["sigma"] <- -T_ / phi["sigma"] - 
      .5 * sum(R_out$variation_out[, "Sigma_diag_grad_sigma"]) / phi["sigma"] ^ 2 +
      (sum(z_vars) + sum_mus_sqrd) / phi["sigma"] ^ 3
    
    ## update gradient by subtracting sum of gradients of E[h(V_t)]
    # for vartheta
    grad_phi["vartheta"] <- grad_phi["vartheta"] -
      sum(E_hV_contrib$integrals_m * 2 * delta ^ (T_:1) * scaled_sums_store_2 / ((1 + delta) * s2 - 1))
    # for nu1, nu2
    for (curVar in c("nu1", "nu2")){
      grad_phi[curVar] <- grad_phi[curVar] -
        sum(
          E_hV_contrib$integrals_m *
            vapply(1:T_, function(ell) scaled_sum(rate = delta, vec = alphasOut$alphas_grad_mat[, curVar][1:ell]), 0)
        )
    }
    # for sigma, tau, beta1, beta2 also add gradients of one half log det Sigma
    for (curVar in c("sigma", "tau", "beta1", "beta2")){
      grad_phi[curVar] <- grad_phi[curVar] -
        sum(
          E_hV_contrib$integrals_s *
            # this is the partial derivative of s_vec wrt the current variable (rho, sigma or tau_1)
            R_out$variation_out[, paste0("s_vec_grad_", curVar)]
        ) +
          #partial derivative of log det Sigma
        R_out$det_out[[curVar]]
    }
    
    ## add gradient contributions of E[p(y[T] = 1 | V_T)] = kappa * E_q v_T
    # for vartheta
    grad_phi["vartheta"] <- grad_phi["vartheta"] +
      (K / (1 + delta)) * 2 * delta * (1 + 1 / ((1 + delta) * s2 - 1))
    # for nu1 and nu2 variables
    for (curVar in c("nu1", "nu2")){
      grad_phi[curVar] <- grad_phi[curVar] +
        (K / (1 + delta)) * alphasOut$alphas_grad_mat[, curVar][T_]
    }
    
    retList$grad_phi <- grad_phi
    
    ### NEXT compute gradient wrt W
    grad_W <-
      # start with gradient contribution from E_{q(z)} sum_t log p(z_t)
      ( 1 / (phi["sigma"] ^ 2) ) *
        (scaled_cum_spike_mat[T_, ] + scaled_cum_spike_mat[T_ - 1, ]) * delta * scaled_cumsums_mus[T_] / ((1 + delta) * s2 - 1) -
      # subtract sum of gradients of E[h(V_t)]
      t(scaled_cum_spike_mat -
                matrix(scaled_cum_spike_mat[T_, ] + scaled_cum_spike_mat[T_ - 1, ], ncol = d, nrow = T_, byrow = T) *
                matrix(delta ^ (T_:1) * scaled_sums_store_2, nrow = T_, ncol = d, byrow = F) / ((1 + delta) * s2 - 1)) %*%
        E_hV_contrib$integrals_m +
      # add gradient contribution of E[p(y[T] = 1 | V_T)] = kappa * E_q v_T
      (K / (1 + delta)) * (X[T_ - 1, ] - (scaled_cum_spike_mat[T_, ] + scaled_cum_spike_mat[T_ - 1, ]) * delta / ((1 + delta) * s2 - 1))

    #retList$g_chk <- ( 1 / (phi["sigma"] ^ 2) ) *
    #  (scaled_cum_spike_mat[T_, ] + scaled_cum_spike_mat[T_ - 1, ]) * delta * scaled_cumsums_mus[T_] / ((1 + delta) * s2 - 1)
    
    retList$grad_W <- grad_W
    
    ### LASTLY compute gradient wrt eta
    grad_eta <-
      # add gradient contribution from E_{q(z)} sum_t log p(z_t)
      ( 1 / (phi["sigma"] ^ 2) ) *
        sum(scaled_sums_store[(T_ - 1):T_]) * delta * scaled_cumsums_mus[T_] / ((1 + delta) * s2 - 1) -
      # subtract sum of gradients of E[h(V_t)]
      sum(
        E_hV_contrib$integrals_m *
          (scaled_sums_store - sum(scaled_sums_store[(T_ - 1):T_]) * delta ^ (T_:1) * scaled_sums_store_2 / ((1 + delta) * s2 - 1))
      ) +
      # add gradient contribution of E[p(y[T] = 1 | V_T)] = kappa * E_q v_T
      (K / (1 + delta)) * (1 - sum(scaled_sums_store[(T_ - 1):T_]) * delta / ((1 + delta) * s2 - 1))
    
    retList$grad_eta <- unname(grad_eta)
  }
  
  return(retList)
}

# calculate EM objective and gradient
# over set of indexed spike intervals
# spike intervals indexed by (i, k) for the ith node and its kth observed interspike interval
# function takes cluster object "cl" from parallel package
# and parallelizes the obj_grad_single_spike_out computations over the set of spike intervals
par_obj_grad_spikes <- function(
  cl, interval_indices,
  input_list,
  spike_intervals_df, binary_spike_train_mat,
  delta, K=50, piecewise_approx_list,
  grads = TRUE,
  restrict_W = FALSE,
  W_lb = NULL, W_ub = NULL,
  transform_eta = TRUE,
  eta_lb = NULL, eta_ub = NULL,
  ratio_bound_for_reparam = NULL,
  retListObjs = FALSE,
  sort_on_T = FALSE,
  par_fn = lapply_parLapply
){
  n_interval <- length(interval_indices)
  
  # get permutation of interval indices
  # based on descending order of interval length
  # (idea is we want to run longer intervals first in our parallelized processing)
  if(sort_on_T){
    interval_indices <-
      spike_intervals_df %>%
      filter(index %in% interval_indices) %>%
      arrange(desc(T_)) %>%
      pull(index)
  }
  
  obj_grad_list <- par_fn(
    cl,
    interval_indices,
    function(interval_index){
      print(interval_index)
      interval_data <- get_interval_data(interval_index, spike_intervals_df, binary_spike_train_mat)
      
      ## get relevant parameters for node
      W_i <- input_list$W[, interval_data$i]
      eta_i <- input_list$eta[interval_data$i]
      phi_i <- input_list$phi[, interval_data$i]
      
      ## get obj and grad for single spike
      obj_grad_single_spike_out <- obj_grad_single_spike(
        W_in=W_i, phi=phi_i, eta=eta_i,
        T_=interval_data$T_, X=interval_data$X,
        delta = delta,
        K=50, piecewise_approx_list=piecewise_approx_list,
        grads = grads
      )
      
      if (grads){
        retList <- c(
          list(obj_vec = numeric(ncol(input_list$W))),
          lapply(
            input_list, function(input_var) input_var * 0
          )
        )
        retList$obj_vec[interval_data$i] <- obj_grad_single_spike_out$obj
        
        retList$phi[, interval_data$i] <- obj_grad_single_spike_out$grad_phi
        retList$eta[interval_data$i] <- obj_grad_single_spike_out$grad_eta
        retList$W[, interval_data$i] <- obj_grad_single_spike_out$grad_W
      } else {
        retList <- list(
          obj = obj_grad_single_spike_out$obj
        )
      }
      
      return(retList)
    }
  )
  
  if (retListObjs){
    if (grads){
      return(obj_grad_list)
    } else {
      return(unlist(obj_grad_list))
    }
  }
  
  summed_list <- list()
  for (cur_name in names(obj_grad_list[[1]])){
    summed_list[[cur_name]] <- Reduce("+", lapply(obj_grad_list, "[[", cur_name))
  }
  
  
  # get gradients wrt to the reparameterized phi variables
  if (grads) {
    names(summed_list)[2:4] <- paste0("grad_", names(summed_list)[2:4])
    
    # phi vars
    phi_grad_reparam <- sapply(
      1:ncol(input_list$W),
      function(node_i){
        get_grad_reparam_phi(
          grad_phi = summed_list$grad_phi[, node_i],
          phi = input_list$phi[, node_i],
          ratio_bound_for_reparam = ratio_bound_for_reparam
        )
      }
    )
    summed_list$grad_phi <- phi_grad_reparam
    # gradient reparameterization for logit tranformation of "W" params if bounds imposed
    if (restrict_W){
      summed_list$grad_W <- summed_list$grad_W * (input_list$W - W_lb) * (W_ub - input_list$W) / (W_ub - W_lb)
    }
    # and for eta if bounds imposed
    if(transform_eta){
      summed_list$grad_eta <- summed_list$grad_eta * (input_list$eta - eta_lb) * (eta_ub - input_list$eta) / (eta_ub - eta_lb)
    }
  }
  
  for (cur_name in names(summed_list)){
    summed_list[[paste0(cur_name, "_", "mean")]] <- summed_list[[cur_name]] / n_interval
    if (!(cur_name %in% c("obj", "obj_vec"))){
      summed_list[[paste0("sum_sqrd_", cur_name, "_", "mean")]] <- sum(summed_list[[cur_name]] ^ 2) / n_interval
    }
  }
  
  return(summed_list)
}

# function for performing adam update steps of ALL variables
update_adam_all <- function(
  input_list,
  training_indices, #the interval indices used to update phi
  spike_intervals_df, binary_spike_train_mat,
  delta, K=50, piecewise_approx_list,
  numNodesPar,
  minibatch_count, # set minbatch sizes in function
  alpha_step_size,
  m_v_lists_init,
  num_epochs_run=1,
  beta_1=0.9, beta_2=0.999, epsilon=1e-8,
  restrict_W = FALSE,
  W_lb = -.2, W_ub = .2,
  eta_lb=-.1, eta_ub=.2, transform_eta=TRUE,
  ratio_bound_for_reparam=NULL,
  restict_print_cols=NULL,
  epochs_save = FALSE,
  T_train_cutoff = NULL
){
  if(epochs_save){
    if (is.null(T_train_cutoff)){
      epochs_save_suffix <- stringr::str_extract_all(Sys.time(), "[0-9]") %>% unlist() %>% paste(collapse = '')
    } else {
      epochs_save_suffix <-
        paste0(T_train_cutoff, "_", stringr::str_extract_all(Sys.time(), "[0-9]") %>% unlist() %>% paste(collapse = ''))
    }
  }
  
  n <- ncol(input_list$W)
  num_intervals <- length(training_indices)
  base_minibatch_size = floor(num_intervals / minibatch_count)
  remainder <- num_intervals - minibatch_count * base_minibatch_size
  minibatch_sizes = c(rep(base_minibatch_size, minibatch_count - remainder), rep(base_minibatch_size + 1, remainder))
  m_v_lists <- m_v_lists_init
  
  param_list = list(
    ratio_bound_for_reparam=ratio_bound_for_reparam,
    eta_lb=eta_lb,
    eta_ub=eta_ub,
    restrict_W = restrict_W,
    W_lb = W_lb,
    W_ub = W_ub
  )
  
  # keep track of phis, gradients and objectives at epochs
  values_epoch_store <- vector(mode = "list", length = num_epochs_run + 1)
  values_epoch_store[[1]] <- input_list
  epoch_objs_mat <- matrix(NA, nrow = length(train_indices), ncol = num_epochs_run + 1)
  obj_epoch_train_store <- numeric(num_epochs_run + 1)
  epoch_runtimes_store <- numeric(num_epochs_run)
  
  cl <- makeCluster(numNodesPar, type="FORK")
  cur_epoch <- 0
  
  startTime <- Sys.time()
  cat("initial obj: ")
  epoch_objs_mat[, cur_epoch + 1] <- par_obj_grad_spikes(
    cl=cl,
    input_list=input_list,
    interval_indices=training_indices,
    spike_intervals_df=spike_intervals_df, binary_spike_train_mat=binary_spike_train_mat,
    delta=delta, K=K, piecewise_approx_list=piecewise_approx_list,
    grads = FALSE, retListObjs = TRUE
  )
  obj_epoch_train_store[cur_epoch + 1] <- sum(epoch_objs_mat[, cur_epoch + 1])
  cat(obj_epoch_train_store[cur_epoch + 1], "obj runtime: ")
  print(Sys.time() - startTime)
  
  while(cur_epoch < num_epochs_run){
    epochStoreTime <- Sys.time()
    cur_epoch <- cur_epoch + 1
    print(paste("epoch", cur_epoch, "of", num_epochs_run))
    
    # shuffle intervals for random minibatch assignments
    intervals_shuffle <- sample(training_indices, length(training_indices))
    
    for (minibatch_index in 1:minibatch_count){
      # print updates
      if (TRUE){
        #print(Sys.time() - epochStoreTime)
        #print(paste("running update on minibatch", minibatch_index, "of", minibatch_count,
        #            "in epoch", cur_epoch, "of", num_epochs_run))
        if (minibatch_index == 1){
          if(is.null(restict_print_cols)){
            print(rbind(input_list$phi, input_list$eta))
            cat("running", minibatch_count, "minibatch updates in epoch", cur_epoch, "of", num_epochs_run, "\n")
          } else {
            print(rbind(input_list$phi[, restict_print_cols], input_list$eta[restict_print_cols]))
          }
        }
        cat(minibatch_index, "")
      }

      minibatch_shuffled_interval_indices <-
        (minibatch_index > 1) * sum(minibatch_sizes[1:(minibatch_index-1)]) +
        (1:minibatch_sizes[minibatch_index])
      minibatch_interval_indices <- intervals_shuffle[minibatch_shuffled_interval_indices]
      
      # run main function to calculate objective and gradients with current phi
      obj_grad_minibatch <- par_obj_grad_spikes(
        cl=cl,
        input_list=input_list,
        interval_indices=minibatch_interval_indices,
        spike_intervals_df=spike_intervals_df, binary_spike_train_mat=binary_spike_train_mat,
        delta=delta, K=K, piecewise_approx_list=piecewise_approx_list,
        grads = TRUE,
        restrict_W=restrict_W, W_lb=W_lb, W_ub=W_ub,
        eta_lb=eta_lb, eta_ub=eta_ub, transform_eta=TRUE,
        sort_on_T = FALSE,
        par_fn = clusterApplyLB
      )
      
      # adam updates
      # using mean gradients
      for (var_name in names(input_list)){
        m <-
          beta_1 * m_v_lists$m[[var_name]] + (1 - beta_1) * obj_grad_minibatch[[paste0("grad_", var_name, "_mean")]]
        v <-
          beta_2 * m_v_lists$v[[var_name]] + (1 - beta_2) * obj_grad_minibatch[[paste0("grad_", var_name, "_mean")]] ^ 2
        
        # store for future use
        m_v_lists$m[[var_name]] <- m
        m_v_lists$v[[var_name]] <- v
        
        # get bias corrected bits for udate
        m_bias_corrected <- m / (1 - beta_1)
        v_bias_corrected <- v / (1 - beta_2)
        
        var_reparam_before_update <-
          get_reparam_from_variable(input_list[[var_name]], var_name, param_list)
        
        var_reparam <- var_reparam_before_update +
          alpha_step_size * m_bias_corrected / (sqrt(v_bias_corrected) + epsilon)
        
        var_updated <-
          get_variable_from_reparam(var_reparam, var_name, param_list)
        
        input_list[[var_name]] <- var_updated
      }
      
      clusterExport(cl, "input_list")
    }
    
    values_epoch_store[[cur_epoch + 1]] <- input_list
    
    cat("\ncalc obj cur epoch: ")
    #print(Sys.time() - epochStoreTime)
    epoch_objs_mat[, cur_epoch + 1] <- par_obj_grad_spikes(
      cl=cl,
      input_list=input_list,
      interval_indices=training_indices,
      spike_intervals_df=spike_intervals_df, binary_spike_train_mat=binary_spike_train_mat,
      delta=delta, K=K, piecewise_approx_list=piecewise_approx_list,
      grads = FALSE, retListObjs = TRUE
    )
    obj_epoch_train_store[cur_epoch + 1] <- sum(epoch_objs_mat[, cur_epoch + 1])
    cat(obj_epoch_train_store[cur_epoch + 1])
    cat(" epoch runtime:", Sys.time() - epochStoreTime, "\n")
    epoch_runtimes_store[cur_epoch] <- Sys.time() - epochStoreTime
    
    if (epochs_save) {
      saveRDS(
        list(
          epoch_runtimes_store = epoch_runtimes_store,
          values_epoch = values_epoch_store,
          epoch_objs_mat = epoch_objs_mat,
          obj_epoch_train = obj_epoch_train_store,
          m_v_lists_final = m_v_lists,
          values_out = input_list,
          epochs_completed = cur_epoch
        ),
        file=paste0("xx - output/epochs_save/epochs_save_", epochs_save_suffix, ".RData")
      )
    }
  }
  stopCluster(cl)
  
  return_list <-
    list(
      epoch_runtimes_store = epoch_runtimes_store,
      values_epoch = values_epoch_store,
      epoch_objs_mat = epoch_objs_mat,
      obj_epoch_train = obj_epoch_train_store,
      m_v_lists_final = m_v_lists,
      values_out = input_list,
      epochs = num_epochs_run
    )
  
  return(return_list)
}
