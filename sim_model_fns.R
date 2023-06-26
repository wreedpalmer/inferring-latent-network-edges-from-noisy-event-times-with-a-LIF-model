# sample multivariate Gaussian z given mean vector mus
# along with diagonal (Sigma11) and off diagonal (Sigma12) entries of double constant covariance matrix
# inputs:
#   T_ - spike interval length
#   a list fitVals with the following named components:
#     mus - mean vectoir
#     Sigma11 - diagonal element of covariance
#     Sigma12 - off diagonal element of covariance
#     non_z_cum_input_vec - for constructing voltage path
#   delta - leakage parameter
#   return_z - boolean whether to return sampled z's or voltage path
sample_from_fit <- function(fitVals, delta, return_z = TRUE){
  T_ <- length(fitVals$mus)
  z_out <- numeric(T_)
  
  # iid std normal sample for construction of zs
  epsilon_vec <- rnorm(T_)
  
  # vector of the conditional variance of z_t | z_1,...,z_r for t > r
  c_var_vec <- numeric(T_)
  # vector of the conditional covariance of z_s and z_t given z_1,...,z_r for t > s > r
  c_cov_vec <- numeric(T_)
  
  c_var_vec[1] <- fitVals$Sigma_11
  c_cov_vec[1] <- fitVals$Sigma_12
  
  mus <- fitVals$mus
  mus_use <- numeric(T_)
  mus_use[1] <- mus[1]
  
  z_out[1] <- mus[1] + epsilon_vec[1] * sqrt(c_var_vec[1])
  
  for(k in 2:T_){
    mus <- mus + (z_out[k-1] - mus[k-1]) * c_cov_vec[k - 1] / c_var_vec[k - 1]
    mus_use[k] <- mus[k]
    
    #conditional variance / covariance given z_1,...,z_{k-1}
    c_var_vec[k] <- c_var_vec[k - 1] - c_cov_vec[k - 1] ^ 2 / c_var_vec[k - 1]
    c_cov_vec[k] <- c_cov_vec[k - 1] - c_cov_vec[k - 1] ^ 2 / c_var_vec[k - 1]
    
    z_out[k] <- mus[k] + epsilon_vec[k] * sqrt(c_var_vec[k])
  }
  
  if(return_z){
    return(z_out)
  }
  
  z_scaled_sums <- vapply(1:T_, function(ell) scaled_sum(rate = delta, vec = z_out[1:ell]), 0)
  v_out <- fitVals$non_z_cum_input_vec + z_scaled_sums
  
  return(v_out)
}


##### SIMULATE ON NETWORK
# T_total - total time periods to run
# W - weighted adjacency matrix
# eta_vec - constant input for each node
# sigma_vec - variation in z
# delta - rate of decay
# K - kappa for soft threshold triggering mechanism
runSimNetwork <- function(T_total, warmup,
                          W, V_init=NULL,
                          sigma_vec, eta_vec,
                          delta, K=50){
  
  n <- ncol(W)
  
  outMatList <- list()
  for (i in 1:n){
    outMatList[[i]] <- matrix(0, nrow = T_total - warmup, ncol=11)
    colnames(outMatList[[i]]) <- c("t", "spike_count", "t_ws", "z", "xb", "I_total", "V", "y",
                                   "E_V_w_network", "E_V_no_network", "intercept")
  }
  
  networkPredictors <- list()
  for (i in 1:n){
    networkPredictors[[i]] <- matrix(0, nrow = T_total - warmup, ncol=n)
  }
  
  binarySpikeMat <- matrix(NA, ncol=n, nrow=T_total - warmup)
  
  spikes_mat <- numeric(0)
  spikes_list <- lapply(1:n, function(i) numeric(0))
  
  if (is.null(V_init)){
    V <- numeric(n)
  }else{
    V <- V_init
  }
  
  E_V_w_network <- numeric(n)
  E_V_no_network <- numeric(n)
  intercept <- numeric(n)
  
  net_pred_cur <- matrix(0, nrow = n, ncol = n)
  
  spike_count <- 1
  t_ws <- numeric(n)
  cur_out <- numeric(n)
  
  for (t in 1:T_total){
    #update within spike time count
    t_ws <- t_ws + 1
    
    #get fluctuation
    z <- rnorm(n, 0, 1) * sigma_vec
    
    ## update state
    # include (1 - cur_out) term so that nodes that just spiked do not
    # receive any charge from the most recent set of spikes
    V <- delta * V + eta_vec + z + (cur_out %*% W) * (1 - cur_out)
    
    net_pred_cur <- delta * net_pred_cur + cur_out
    
    E_V_w_network <- delta * E_V_w_network + eta_vec + (cur_out %*% W) * (1 - cur_out)
    E_V_no_network <- delta * E_V_no_network + eta_vec
    intercept <- delta * intercept + 1
    
    #store new row for just 1st node
    #df <- df %>% add_row(t=t, spike_count=spike_count, t_ws=t_ws, xb=sum(betas * x_in), V=V, y=y)
    if(t > warmup){
      for (i in 1:n){
        outMatList[[i]][t - warmup,] <- c(
          t - warmup,
          length(spikes_list[[i]]),
          t_ws[i], z[i],
          ((cur_out %*% W) * (1 - cur_out))[i],              # input from spikes
          eta_vec[i] + ((cur_out %*% W) * (1 - cur_out))[i],        # total input
          V[i], 0,
          E_V_w_network[i], E_V_no_network[i],
          intercept[i]
        )
        
        networkPredictors[[i]][t - warmup,] <-
          net_pred_cur[,i]
      }
    }
    
    # spike or not
    cur_out <- vapply(V, function(v) rbinom(1, 1, gtools::inv.logit(K * (v - 1))), 1)
    
    # store full output
    if(t > warmup){
      binarySpikeMat[t - warmup, ] <- cur_out
      for (i in 1:n){
        outMatList[[i]][t - warmup, "y"] <- cur_out[i]
      }
    }
    
    # reset voltages and within-spike counts for any spikes
    for(i in which(cur_out==1)){
      
      # record spikes after warmup period
      if(t > warmup){
        spikes_list[[i]] <- c(spikes_list[[i]], t)
        spikes_mat <- rbind(spikes_mat, c(spike_count, t, i, t_ws[i],
                                          length(spikes_list[[i]])))
        spike_count <- spike_count + 1
      }
      
      t_ws[i] <- 0
      V[i] <- 0
      E_V_w_network[i] <- 0
      E_V_no_network[i] <- 0
      intercept[i] <- 0
      net_pred_cur[,i] <- 0
    }
  }
  
  #some post processing
  colnames(spikes_mat) <- c("index", "t", "i", "T_", "k")
  
  spike_intervals <- spikes_mat %>%
    as_tibble() %>%
    filter(k > 1) %>%
    mutate(
      t_start = t - T_ + 1 - warmup,
      t_end = t - warmup,
      index = row_number()
    )
  
  return(
    list(
      spikes_list=spikes_list,
      spikes_mat=spikes_mat,
      binary_spike_train=binarySpikeMat,
      outMatList=outMatList,
      spike_intervals=spike_intervals,
      networkPredictors=networkPredictors
    )
  )
}
