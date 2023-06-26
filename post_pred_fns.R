# version of parSapply that will revert to regular lapply if cl is NULL
lapply_parSapply <- function(cl, X, FUN){
  if(is.null(cl)){
    sapply(X, FUN)
  }else{
    parSapply(cl, X, FUN)
  }
}

get_source_index_in_W <- function(source_i, target_j){
  if (source_i < target_j) return(source_i) else return(source_i - 1)
}

# function for getting only the piece of the objective we are interested in
get_component_1 <- function(
  interval_data,
  estimated_parameters,
  delta,
  kappa,
  R_out,
  piecewise_approx_list,
  non_z_cum_input_vec = NULL
){
  T_ <- interval_data$T_
  phi <- estimated_parameters$phi[, interval_data$i]
  knots <- lapply(piecewise_approx_list, "[[", "interval") %>% unlist %>% unique
  
  if(is.null(non_z_cum_input_vec)){
    # cumulative/scaled sum of input *moderated* by W_in
    scaled_cum_arrivals <- c(0, scaled_sum_mat_mult(delta, interval_data$X, estimated_parameters$W[, interval_data$i]))
    
    # total nonrandom input
    non_z_cum_input_vec <- scaled_cum_arrivals +
      estimated_parameters$eta[interval_data$i] * (1 - delta ^ (1:T_)) / (1 - delta)
  }
  
  # get alphas vector
  alphas_vec <- get_alphas(T_, phi["nu1"], phi["nu2"], delta, retGrad = FALSE)$alphas
  
  # save key scaled sum
  s2 <- scaled_sum(rate = delta ^ 2, len = T_)
  
  # expectations E z_t, t=1 to t=T_ under q_phi
  c_mu <- (2 * phi["vartheta"] - non_z_cum_input_vec[T_] - non_z_cum_input_vec[T_ - 1]) /
    ((1 + delta) * s2 - 1)
  mus <- c_mu * delta ^ (T_:1) + alphas_vec
  
  # calculate mean vec E v_t, t=1 to t=T_ under q_phi
  scaled_cumsums_mus <- vapply(1:T_, function(ell) scaled_sum(rate = delta, vec = mus[1:ell]), 0)
  m_vec <- non_z_cum_input_vec + scaled_cumsums_mus
  
  # get s_vec from stored R_outs_list
  #s_vec <- R_outs_list[[interval_data$i]][[paste(T_)]]$variation_out[, "s_vec"]
  s_vec <- R_out$variation_out[, "s_vec"]
  
  R1_mat <- get_R_1_mat_6(mu=m_vec, sigma=s_vec)
  R2_mat_list <- lapply(knots, function(knot) get_R_2_mat_6(x=knot, mu=m_vec, sigma=s_vec))
  pn_list <- lapply(knots, function(knot) pnorm(knot, m_vec, s_vec))
  dn_list <- lapply(knots, function(knot) dnorm(knot, m_vec, s_vec))
  
  neg_obj_contribs <- vapply(
    1:length(piecewise_approx_list),
    function(segment_index){
      segment <- piecewise_approx_list[[segment_index]]
      
      # much more efficient computation
      # see details...
      E_xk_vec_mat <- t(2 * (pn_list[[segment_index + 1]] - pn_list[[segment_index]]) * R1_mat -
                          dn_list[[segment_index + 1]] * R2_mat_list[[segment_index + 1]] +
                          dn_list[[segment_index]] * R2_mat_list[[segment_index]])
      
      # objecitve contribution (negative)
      return(sum(E_xk_vec_mat[1:length(segment$coefs), ] * segment$coefs))
    },
    1
  )
  
  obj_component1 <-
    kappa * (mus[T_] + sum(estimated_parameters$W[, interval_data$i] * interval_data$X[T_ - 1, ])) / (1 + delta) -
    sum(mus ^ 2) / (2 * phi["sigma"] ^ 2) - sum(neg_obj_contribs)
  
  return(
    obj_component1
  )
}

get_component_1_multiple_cols <- function(
    phi,
    delta,
    kappa,
    s_vec_in,
    piecewise_approx_list,
    non_z_cum_input_mat,
    mod_input_2nd_last_period
){
  T_ <- length(s_vec_in)
  knots <- lapply(piecewise_approx_list, "[[", "interval") %>% unlist %>% unique
  
  # save key scaled sum
  s2 <- scaled_sum(rate = delta ^ 2, len = T_)
  
  # expectations E z_t, t=1 to t=T_ under q_phi
  c_mu_vec <- (2 * phi["vartheta"] - non_z_cum_input_mat[T_, ] - non_z_cum_input_mat[T_ - 1, ]) /
    ((1 + delta) * s2 - 1)
  mu_mat <- get_alphas(T_, phi["nu1"], phi["nu2"], delta, retGrad = FALSE)$alphas + delta ^ (T_:1) %*% t(c_mu_vec)
  
  # calculate mean mat E v_t, t=1 to t=T_ under q_phi
  scaled_cumsums_mu_mat <- scaled_sum_mat(delta, mu_mat)
  m_mat <- non_z_cum_input_mat + scaled_cumsums_mu_mat
  m_mat_as_vec <- as.numeric(m_mat)
  
  s_vec <- rep(s_vec_in, ncol(non_z_cum_input_mat))
  
  R1_mat <- get_R_1_mat_6(mu=m_mat_as_vec, sigma=s_vec)
  R2_mat_list <- lapply(knots, function(knot) get_R_2_mat_6(x=knot, mu=m_mat_as_vec, sigma=s_vec))
  pn_list <- lapply(knots, function(knot) pnorm(knot, m_mat_as_vec, s_vec))
  dn_list <- lapply(knots, function(knot) dnorm(knot, m_mat_as_vec, s_vec))
  
  neg_obj_contribs_unsummed <- sapply(
    1:length(piecewise_approx_list),
    function(segment_index){
      segment <- piecewise_approx_list[[segment_index]]
      
      # much more efficient computation
      # see details...
      E_xk_vec_mat <- t(2 * (pn_list[[segment_index + 1]] - pn_list[[segment_index]]) * R1_mat -
                          dn_list[[segment_index + 1]] * R2_mat_list[[segment_index + 1]] +
                          dn_list[[segment_index]] * R2_mat_list[[segment_index]])
      
      # objecitve contribution (negative)
      return(colSums(E_xk_vec_mat[1:length(segment$coefs), ] * segment$coefs))
    }
  ) %>%
    rowSums()
  
  neg_obj_contribs_vec <- colSums(matrix(neg_obj_contribs_unsummed, ncol = ncol(non_z_cum_input_mat)))
  
  obj_component_vec <-
    kappa * (mu_mat[T_, ] + mod_input_2nd_last_period) / (1 + delta) -
    colSums(mu_mat ^ 2) / (2 * phi["sigma"] ^ 2) - neg_obj_contribs_vec
  
  return(
    setNames(obj_component_vec, colnames(non_z_cum_input_mat))
  )
}


get_s_vec_mats_list <- function(
    runList,
    T_train_cutoff,
    T_test_cutoff,
    sim = TRUE,
    parallelize_over_nodes = TRUE,
    n_nodes_par = NULL
){
  if (sim) {
    spike_intervals_df <- runList$model_params$sim$spike_intervals
    delta <- runList$model_params$delta
  } else {
    spike_intervals_df <- runList$model_params$spike_intervals_df %>%
      filter(filter_length == "keep")
    delta <- runList$fit_params$delta
  }
  
  nodes <- 1:runList$model_params$n
  
  last_observed_spike_times <- spike_intervals_df %>%
    filter(t_end <= T_train_cutoff & i %in% nodes) %>%
    group_by(i) %>%
    summarize(last_spike = max(t_end))
  
  node_length_pairs <- spike_intervals_df %>%
    filter(
      t_end > min(last_observed_spike_times$last_spike) + 1 &
        t_start < T_test_cutoff &
        i %in% nodes) %>%
    group_by(i, T_) %>%
    summarize(count = n())
  
  max_Ts <-
    node_length_pairs %>%
    group_by(i) %>%
    summarize(max_T = max(T_))
  
  if(parallelize_over_nodes){
    cl <- makeCluster(n_nodes_par, type="FORK")
  } else {
    cl <- NULL
  }
  s_vec_mats <-
    lapply_parLapplyLB(
      cl,
      nodes,
      function(cur_node){
        print(cur_node)
        phi <- runList$adam_update$values_out$phi[, cur_node]
        T_s <- 1:max_Ts$max_T[cur_node]
        out_row <- numeric(max_Ts$max_T[cur_node])
        s_vec_mat <- sapply(
          T_s, function(T_){
            #print(T_)
            c_vec <- .5 * (delta ^ ((T_ - 1):0) + c(delta ^ ((T_ - 2):0), 0))
            s_vec <-
              get_R(
                T_ = T_,
                phi = phi,
                delta = delta,
                c_vec = c_vec,
                grads = F,
                get_s_vec_out = TRUE
              )$variation_out[,"s_vec"]
            out_row[1:T_] <- s_vec
            return(out_row)
          }
        ) %>% t
        return(s_vec_mat)
      }
    )
  names(s_vec_mats) <- paste(nodes)
  
  if(parallelize_over_nodes){
    stopCluster(cl)
  }
  
  return(s_vec_mats)
}

######## get pre-stored R_outs
store_R_outs <- function(
  runList,
  T_train_cutoff,
  T_test_cutoff,
  cl = NULL,
  sim = TRUE,
  nodesList = NULL,
  listOflists = TRUE
){
  if (is.null(runList$model_params$delta)){
    delta <- runList$fit_params$delta
  } else {
    delta <- runList$model_params$delta
  }
  
  if(!is.null(cl)){
    clusterExport(cl,"delta",environment())
  }
  
  get_R_lapply <- function(T_){
    c_vec <- .5 * (delta ^ ((T_ - 1):0) + c(delta ^ ((T_ - 2):0), 0))
    R_out <- get_R(T_ = T_, phi = phi, delta = delta, c_vec = c_vec, grads = F, get_s_vec_out = TRUE)
    return(R_out)
  }
  environment(get_R_lapply) <- .GlobalEnv
  
  n <- runList$model_params$n
  if(is.null(nodesList)){
    nodesList <- 1:n
  }
  
  if (sim) {
    spike_intervals_df <- runList$model_params$sim$spike_intervals
  } else {
    spike_intervals_df <- runList$model_params$spike_intervals_df %>% filter(filter_length == "keep")
  }
  
  last_observed_spike_times <- spike_intervals_df %>%
    filter(t_end <= T_train_cutoff & i %in% nodesList) %>%
    group_by(i) %>%
    summarize(last_spike = max(t_end))
  
  node_length_pairs <- spike_intervals_df %>%
    filter(
      t_end > min(last_observed_spike_times$last_spike) + 1 &
        t_start < T_test_cutoff &
        #t_end <= T_test_cutoff &
        i %in% nodesList) %>%
    group_by(i, T_) %>%
    summarize(count = n())
  
  if (listOflists) {
    R_outs_store <-
      lapply(nodesList, function(cur_node){
        print(cur_node)
        T_s <- node_length_pairs %>% filter(i == cur_node) %>% pull(T_)
        phi <- runList$adam_update$values_out$phi[, cur_node]
        clusterExport(cl,"phi",environment())
        out <- setNames(
          lapply_parLapply(
            cl=cl,
            T_s, function(T_){
              print(T_)
              get_R_lapply(T_)
            }
          ),
          as.character(T_s)
        )
        return(out)
      })
    names(R_outs_store) <- paste(nodesList)
    
    return(R_outs_store)
  } else{
    max_Ts <- #spike_intervals_df %>%
      node_length_pairs %>%
      group_by(i) %>%
      summarize(max_T = max(T_))
    
    s_vec_mats <-
      lapply_parLapplyLB(
        cl,
        nodesList,
        function(cur_node){
          #print(cur_node)
          phi <- runList$adam_update$values_out$phi[, cur_node]
          T_s <- 1:max_Ts$max_T[cur_node]
          out_row <- numeric(max_Ts$max_T[cur_node])
          s_vec_mat <- sapply(
            T_s, function(T_){
              #print(T_)
              c_vec <- .5 * (delta ^ ((T_ - 1):0) + c(delta ^ ((T_ - 2):0), 0))
              s_vec <-
                get_R(
                  T_ = T_,
                  phi = phi,
                  delta = delta,
                  c_vec = c_vec,
                  grads = F,
                  get_s_vec_out = TRUE
                )$variation_out[,"s_vec"]
              out_row[1:T_] <- s_vec
              return(out_row)
            }
          ) %>% t
          return(s_vec_mat)
        }
      )
    names(s_vec_mats) <- paste(nodesList)
    return(s_vec_mats)
  }
}

# assume voltage is reset reset (spike at prev time)
# draw n_sims (partial or complete) inter-spike intervals, with length at most equal to the length of the provided non_z_cum_input vector
# given parameters and input
# stop each sim as soon as first spike recorded or when max length reached
sim_spike <- function(
  n_sims,
  delta,
  kappa,
  sigma,
  eta,
  non_z_cum_input,
  return_prop_non_spike = FALSE
){
  
  length_out <- length(non_z_cum_input)
  
  first_spikes <- numeric(n_sims)
  no_spike_yet <- rep(TRUE, n_sims)
  rand_comp <- sigma * rnorm(n_sims)
  for (k in 1:length_out){
    if (k > 1){
      rand_comp <- delta * rand_comp + sigma * rnorm(n_sims)
    }
    
    inv_logit_input <- kappa * (rand_comp + non_z_cum_input[k] - 1)
    
    #check_spike <- (runif(n_sims) < exp(inv_logit_input) / (1 + exp(inv_logit_input)))
    
    first_spike <- no_spike_yet & # make sure no spike yet
      (runif(n_sims) < exp(inv_logit_input) / (1 + exp(inv_logit_input)))   # check if spike observed
    first_spikes[first_spike] <- k
    no_spike_yet <- no_spike_yet & !first_spike
    
    if (sum(no_spike_yet) == 0) break
  }
  
  spike_frequencies <- tibble(T_ = first_spikes) %>%
    group_by(T_) %>%
    summarize(freq = n() / n_sims) %>%
    filter(T_ > 0)
  
  ret_vec <- numeric(length_out)
  ret_vec[spike_frequencies$T_] <- spike_frequencies$freq
  
  if (return_prop_non_spike) {
    return(1 - sum(ret_vec))
  } else {
    return(ret_vec)
  }
}

sim_spike_multiple_cols <- function(
    n_sims,
    delta,
    kappa,
    sigma,
    non_z_cum_input
){
  num_input_cols <- ncol(non_z_cum_input)
  length_out <- nrow(non_z_cum_input)
  first_spikes <- matrix(0, nrow = n_sims, ncol = num_input_cols)
  no_spike_yet <- matrix(TRUE, nrow = n_sims, ncol = num_input_cols)
  rand_comp <- sigma * rnorm(n_sims)
  
  for (k in 1:length_out){
    if (k > 1){
      rand_comp <- delta * rand_comp + sigma * rnorm(n_sims)
    }
    # make sure no spike yet & check if spike observed
    first_spike <- no_spike_yet &
      (gtools::logit(runif(n_sims)) <
         sapply(1:num_input_cols, function(col_j) kappa * (rand_comp + non_z_cum_input[k, col_j] - 1)))
    first_spikes[first_spike] <- k
    no_spike_yet <- no_spike_yet & !first_spike
    #if (sum(no_spike_yet) == 0) break
  }
  
  ret_vec <- setNames(numeric(length_out), paste(1:length_out))
  out_mat <-
    sapply(
      1:num_input_cols,
      function(col_j){
        freq <- table(c(0, first_spikes[, col_j]))[-1] / n_sims
        ret_vec[names(freq)] <- freq
        return(unname(ret_vec))
      }
    )
  return(matrix(out_mat, ncol=num_input_cols))
}

# efficiently sample from non-anticipation distribution
sim_no_spikes_yet_SAPPLY <- function(
    n_sims,
    delta,
    kappa,
    sigma,
    non_z_cum_input # T by p matrix
){
  
  fn_environment <- environment()
  num_input_cols <- ncol(non_z_cum_input)
  length_to_check <- nrow(non_z_cum_input)
  no_spike_yet <- matrix(1, nrow = n_sims, ncol = num_input_cols)
  rand_comp <- sigma * rnorm(n_sims)
  ret_row <- rep(1, num_input_cols)
  
  return(
    t(
      sapply(
        1:length_to_check,
        function(t){
          # if all sims and cols have spiked we no longer need to perform any updates
          if (sum(ret_row) > 0){
            if (t > 1){
              ## recursively update scaled random component of latent voltage
              assign("rand_comp", delta * rand_comp + sigma * rnorm(n_sims), envir = fn_environment)
            }
            
            # combine random component with non-random and check probabilistic spiking mechanism
            # update no_spike_yet matrix with observed spikes
            assign("no_spike_yet",
                   no_spike_yet &
                     (gtools::logit(runif(n_sims)) >
                        kappa * (Rfast::Outer(non_z_cum_input[t, ], rand_comp, oper = "+") - 1)),
                   envir = fn_environment)
            
            # return proportion of simulations with no spike yet for each column
            assign("ret_row", colMeans(no_spike_yet), envir = fn_environment)
          }
          return(ret_row)
        }
      )
    )
  )
}

# get simulated spike trains for held out node in testing period
get_NA_sims <- function(
    n_sims,
    X_minus_i,
    delta,
    kappa,
    W_in,
    eta,
    sigma,
    epsilon_fill = 1e-06
){
  
  testing_length <- nrow(X_minus_i)
  
  SIM_MAT_COND <- matrix(0, nrow = testing_length, testing_length)
  SIM_MAT_JOINT <- matrix(0, nrow = testing_length, testing_length)
  PROB_SPIKE_VEC <- numeric(testing_length)
  
  # get non-randum scaled cumulative inputs starting at each time in testing period
  non_z_cum_input_vec <- numeric(testing_length)
  non_z_cum_input_columns <- sapply(
    testing_length:1,
    function(T_){
      if (T_ > 1){
        non_z_cum_input_vec[1:T_] <-
          c(0,
            scaled_sum_mat_mult(
              delta,
              X_minus_i[(testing_length - T_ + 1):(testing_length - 1), , drop = F],
              W_in)
          ) + eta * (1 - delta ^ (1:T_)) / (1 - delta)
      } else {
        non_z_cum_input_vec[1:T_] <- eta * (1 - delta ^ (1:T_)) / (1 - delta)
      }
      return(non_z_cum_input_vec)
    }
  )
  
  # *efficiently* simulate spikes from reset in each time in testing period
  no_spike_yet <- matrix(TRUE, nrow = n_sims, ncol = testing_length)
  rand_comp <- sigma * rnorm(n_sims)
  
  no_spike_yet_MAT <- matrix(0, nrow = testing_length, ncol = testing_length)
  
  for (t in 1:testing_length){
    print(t)
    
    if (t > 1){
      rand_comp <- delta * rand_comp + sigma * rnorm(n_sims)
    }
    
    # check if spike observed
    no_spike_yet <- no_spike_yet[, 1:(testing_length - t + 1), drop = F] &
      (gtools::logit(runif(n_sims)) >
         kappa * (Rfast::Outer(non_z_cum_input_columns[t, 1:(testing_length - t + 1), drop = F], rand_comp, oper = "+") - 1))
    
    ret_row_vals <- colMeans(no_spike_yet)
    
    no_spike_yet_MAT[t, 1:(testing_length - t + 1)] <- ret_row_vals
    
    if (sum(ret_row_vals) == 0) break
  }
  
  # fill in matrices to return
  # perform some probabilistic calculations...
  first_row <- -diff(c(1, no_spike_yet_MAT[, 1]))
  sum_first_row <- sum(first_row)
  first_row[which(first_row == 0 &
                    no_spike_yet_MAT[, 1] > 0) #& 1:testing_length > 1
  ] <- epsilon_fill
  first_row <- first_row * sum_first_row / sum(first_row)
  
  SIM_MAT_COND[1, ] <- first_row -> SIM_MAT_JOINT[1, ]
  #SIM_MAT_COND[1, ] <- -diff(c(1, no_spike_yet_MAT[, 1])) -> SIM_MAT_JOINT[1, ]
  for (start_sim in 2:(testing_length - 1)){
    
    prob_spike_at_prev_time <- sum(SIM_MAT_JOINT[1:(start_sim - 1), start_sim - 1])
    PROB_SPIKE_VEC[start_sim - 1] <- prob_spike_at_prev_time
    
    spike_sim_vec <- -diff(c(1, no_spike_yet_MAT[1:(testing_length - start_sim + 1), start_sim]))
    
    SIM_MAT_JOINT[start_sim, start_sim:testing_length] <- prob_spike_at_prev_time * spike_sim_vec
    SIM_MAT_COND[start_sim, start_sim:testing_length] <- spike_sim_vec
  }
  
  PROB_SPIKE_VEC[(testing_length - 1):testing_length] <- colSums(SIM_MAT_JOINT[, (testing_length - 1):testing_length])
  
  # matrix of conditional probablitities calculated from simulations
  # COND_PROB_1[s, t] = P(y_t = 1 | y_s = 1)
  COND_PROB_1 <- diag(testing_length)
  for(s in 1:testing_length){
    if (PROB_SPIKE_VEC[s] > 0){
      if (s < testing_length){
        #fill in current row above diagonal
        for (t in (s + 1):testing_length){
          COND_PROB_1[s, t] <- sum(SIM_MAT_COND[(s + 1):t, t] * COND_PROB_1[s, s:(t - 1)])
        }
      }
      if (s > 2){
        # fill in next row below diagonal
        COND_PROB_1[s, 1:(s - 1)] <- (COND_PROB_1[1:(s - 1), s] * PROB_SPIKE_VEC[1:(s - 1)]) / PROB_SPIKE_VEC[s]
      }
    }
  }
  
  # matrix of conditional probablitities calculated from simulations
  # COND_PROB_1[s, t] = P(y_t = 1 | y_s = 0)
  # note diag is zero (obviously)
  COND_PROB_0 <- matrix(0, nrow = testing_length, testing_length)
  for(s in 1:(testing_length)){
    if (s < testing_length){
      #fill in current row above diagonal
      for (t in (s + 1):testing_length){
        COND_PROB_0[s, t] <- sum(SIM_MAT_JOINT[1:s, t]) / (1 - PROB_SPIKE_VEC[s]) +
          sum(SIM_MAT_COND[(s + 1):t, t] * COND_PROB_0[s, s:(t - 1)])
      }
    }
    if (s > 2){
      # fill in next row below diagonal
      COND_PROB_0[s, 1:(s - 1)] <- ((1 - COND_PROB_1[1:(s - 1), s]) * PROB_SPIKE_VEC[1:(s - 1)]) / (1 - PROB_SPIKE_VEC[s])
    }
  }
  
  
  return(
    list(
      SIM_MAT_JOINT = SIM_MAT_JOINT,
      SIM_MAT_COND = SIM_MAT_COND,
      PROB_SPIKE_VEC = PROB_SPIKE_VEC,
      COND_PROB_1_MAT = COND_PROB_1,
      COND_PROB_0_MAT = COND_PROB_0
    )
  )
}

run_NA_sims_nodes <- function(
  runList,
  nodes_to_run = NULL,
  n_sims = 10000,
  T_train_cutoff,
  T_test_cutoff,
  parallelize_over_nodes = TRUE,
  sim = TRUE,
  n_nodes_par = NULL
){
  
  if (is.null(nodes_to_run)) {
    nodes_to_run <- 1:runList$model_params$n
  } else {
    nodes_to_run <- intersect(nodes_to_run, 1:runList$model_params$n)
  }
  
  if (sim) {
    spike_intervals_df <- runList$model_params$sim$spike_intervals
    full_spike_mat <- runList$model_params$sim$binary_spike_train
    delta <- runList$model_params$delta
    kappa <- runList$model_params$K
  } else {
    spike_intervals_df <- runList$model_params$spike_intervals_df
    full_spike_mat <- runList$model_params$binary_spike_mat
    delta <- runList$fit_params$delta
    kappa <- runList$fit_params$K
  }
  
  estimated_parameters <- runList$adam_update$values_out
  
  last_observed_spike_times <-
    spike_intervals_df %>%
    filter(t_end <= T_train_cutoff) %>%
    plyr::daply("i", function(df) max(df$t_end))
  
  if(parallelize_over_nodes){
    cl <- makeCluster(n_nodes_par, type="FORK")
  } else {
    cl <- NULL
  }
  outList <-
    lapply_parLapply(
      cl,
      nodes_to_run,
      function(held_out_node){
        retList <- vector(mode = "list", length = 4)
        names(retList) <- c("held_out_node", "last_observed_spike_time", "testing_range", "NA_sims")
        retList$held_out_node <- held_out_node
        
        last_observed_spike_time <- last_observed_spike_times[held_out_node] -> retList$last_observed_spike_time
        testing_range <- (last_observed_spike_time + 1):T_test_cutoff -> retList$testing_range
        
        retList$NA_sims <- get_NA_sims(
          n_sims = n_sims,
          X_minus_i = full_spike_mat[testing_range, -held_out_node],
          delta = delta,
          kappa = kappa,
          W_in = estimated_parameters$W[, held_out_node],
          eta = estimated_parameters$eta[held_out_node],
          sigma = estimated_parameters$phi["sigma", held_out_node]
        )
        
        return(retList)
      }
    )
  if(parallelize_over_nodes){
    stopCluster(cl)
  }
  
  names(outList) <- paste(nodes_to_run)
  outList$n_sims <- n_sims

  outList$spike_intervals_df <- spike_intervals_df
  outList$full_spike_mat <- full_spike_mat
  #outList$last_observed_spike_times <- last_observed_spike_times
  outList$estimated_parameters <- estimated_parameters
  #outList$min_Ts <- spike_intervals_df %>%
  #  group_by(i) %>%
  #  summarize(min_T = min(T_)) %>%
  #  pull(min_T)
  
  outList$T_train_cutoff = T_train_cutoff
  outList$T_test_cutoff = T_test_cutoff
  
  return(outList)
}

# helper function for taking logs and avoiding blowups / neg infinity
log_epsilon_if_small_else_x <- function(x, epsilon){
  ifelse(x < epsilon, log(epsilon), log(x))
}

# helper function for taking logs and avoiding blowups / neg infinity
get_adj_scores <- function(liks_1_2, n_sims, sim_ob_threshold, adj_probs, ret_diff = FALSE){
  obs_1_2 = liks_1_2 * n_sims
  if (max(obs_1_2) < sim_ob_threshold){
    scores_out <- c(0, 0)
  } else {
    scores_out <-
      vapply(1:2, function(ell) ifelse(obs_1_2[ell] < sim_ob_threshold,
                                       log(adj_probs[obs_1_2[ell] + 1]),
                                       log(liks_1_2[ell])),
             0)
  }
  if (ret_diff){
    return(scores_out[1] - scores_out[2])
  } else {
    return(scores_out)
  }
}

# another helper function that we will use twice below
get_adj_score <- function(row, n_sims, sim_ob_threshold, adj_probs, ret_diff = FALSE){
  if(!is.na(row[3])){
    if (min(max(row[1:2]), max(row[3:4])) * n_sims < sim_ob_threshold){
      if(ret_diff) return(0) else return(c(0, 0))
    } else {
      score_adj <- get_adj_scores(row[1:2], n_sims, sim_ob_threshold, adj_probs, ret_diff = FALSE)
      log_p_adj <- get_adj_scores(row[3:4], n_sims, sim_ob_threshold, adj_probs, ret_diff = FALSE)
      ret_row <- pmin(score_adj, log_p_adj) - log_p_adj
      if(ret_diff) return(ret_row[1] - ret_row[2]) else return(ret_row)
    }
  } else {
    return(get_adj_scores(row[1:2], n_sims, sim_ob_threshold, adj_probs, ret_diff = ret_diff))
  }
}

get_fwrd_bckwrd <- function(
    held_out_node,
    sims,
    testing_range,
    estimated_parameters,
    #s_vec_mats_list,
    intervals_testing_excl_held_out,
    binary_spike_mat_testing_period, #binary_spike_mat_testing_period <- full_spike_mat[testing_range, ]
    delta = delta,
    kappa = kappa,
    n_sims = NULL,
    look_ahead_limit_bckwrd = Inf,
    epsilon_0 = 1e-06, # to prevent zero/neg infinity likelihood calculations, should be no bigger than 1 / nsim...
    epsilon_1 = 1e-25, # for final estimates to enable proper entropy loss calculations
    T_train_vec = NULL,
    save_to_file = TRUE,
    filename_w_path = NULL,
    update_filename_w_path = NULL,
    sim_ob_threshold = 10,
    upper_tail_prob_for_adj = .65
){
  startTime <- Sys.time()
  write(paste("started at", startTime), update_filename_w_path)
  S_ <- length(testing_range)
  
  #modes_run_sim_approx <- "sim"
  
  X_testing_held_out_zeroed <- binary_spike_mat_testing_period
  X_testing_held_out_zeroed[, held_out_node] <- 0
  
  #### NEW WAY OF DEALING WITH SMALL PROBABILITIES ####
  # adjusted probabilitities below a threshold
  #upper_tail_prob <- .65
  #threshold <- 10
  adj_probs <-
    vapply(1:sim_ob_threshold,
           function(x) min(zipfR::Rbeta.inv(upper_tail_prob_for_adj, x, n_sims - x + 1), sim_ob_threshold / n_sims), 0)
  
  
  ################################
  # GET INTERVALS #
  ################################
  
  s_vec <- 0:S_
  # list for storing fwrd scored
  s_rs_list <- lapply(
    s_vec,
    function(s){
      if (s == length(testing_range)) {
        return(NULL)
      } else {
        rs <- which(sims$SIM_MAT_COND[s + 1, ] > 0)
        if (length(rs) == 0){
          return(NULL)
        } else {
          s_to_r_NA_probs = sims$SIM_MAT_COND[s + 1, rs]
          return(
            list(
              rs = rs,
              s_to_r_NA = s_to_r_NA_probs,
              frwrd_score_sim = setNames(numeric(length(rs)), paste(rs)),
              p_next_spike_after_cutoff = 1 - sum(s_to_r_NA_probs),
              score_next_spike_after_cutoff = 0
            )
          )
        }
      }
    }
  )
  names(s_rs_list) <- paste(s_vec)
  
  s_r_max_df <-
    tibble(
      s = s_vec,
      r_max = vapply(lapply(s_rs_list, "[[", "rs"), function(vec) ifelse(is.null(vec), NA, max(vec)), 0),
        ### dont elimate any entries
      p_NA = c(1, sims$PROB_SPIKE_VEC), #c(1, pmax(sims$PROB_SPIKE_VEC, .5 / n_sims)),
      held_out_truth = c(1, binary_spike_mat_testing_period[, held_out_node]),
      epsilon_1 = epsilon_1
    )

  print("getting s_intervals_df")
  s_intervals_df <- mapply(
    function(index, t_start_testing_period, t_end_testing_period){
      #print(which(index == intervals_testing_excl_held_out$index))
      ret_df <- s_r_max_df[1:t_end_testing_period, ] %>%
        filter(
          p_NA > 0 
          #& !is.na(r_max)
          ) %>%
        mutate(
          index = index,
          t_start_testing_period = t_start_testing_period,
          t_end_testing_period = t_end_testing_period,
          # could put something like -10 or -20 here instead of 2 since those excluded are likely to be 
          # partial intervals with little prob of spiking before any r
          include_frwrd =
            !is.na(r_max) &
            t_start_testing_period < r_max - 2 & # could expand this to all intervals
                                                  # when p_next_spike_after_cutoff > 0
                                                    # but this would add a lot...
            s < S_,
          partial_wrt_s = t_start_testing_period < s,
          include_bckwrd =
            t_start_testing_period <= (s + look_ahead_limit_bckwrd) &
            s > 0,
          right_censored = t_end_testing_period > S_
          #include_bckwrd_full = include_bckwrd & !partial_wrt_s
        ) %>%
        filter(include_frwrd | include_bckwrd)
      return(ret_df)
    },
    index = intervals_testing_excl_held_out$index,
    t_start_testing_period = intervals_testing_excl_held_out$t_start_testing_period,
    t_end_testing_period = intervals_testing_excl_held_out$t_end_testing_period,
    SIMPLIFY = FALSE
  ) %>%
    bind_rows()
  
  bckwrd_scores_mat_store <-
    matrix(NA, #ncol = 5 + 4 * length(modes_run_sim_approx)
           ncol = 11, nrow = s_intervals_df %>% filter(include_bckwrd) %>% nrow())
  colnames(bckwrd_scores_mat_store) <-
    c("index", "s", "partial", "lik1", "lik0", "s1", "s0",
      "p1", "p0", "score1", "score0")
  
  intervals_include <-
    intervals_testing_excl_held_out %>%
    filter(index %in% unique(s_intervals_df$index))
  
  ################################
  # LOOP OVER INCLUDED INTERVALS #
  ################################
  msg <- paste("running cals over", nrow(intervals_include), "included intervals")
  print(msg)
  write(paste(msg), update_filename_w_path, append = TRUE)
  
  for(k in seq_along(intervals_include$index)){
    print(Sys.time() - startTime)
    storeTime <- Sys.time()
    
    ##################################
    # SET VARIABLES FOR kth INTERVAL #
    ##################################
    
    interval_index <- intervals_include$index[k]
    t_start_rel_testing <- intervals_include[k, "t_start_testing_period", drop = T]
    t_end_rel_testing <- intervals_include[k, "t_end_testing_period", drop = T]
    right_censored = intervals_include[k, "interval_partially_observed", drop = T]
    i <- intervals_include[k, "i", drop = T]
    
    if (right_censored) {
      T_ <- S_ - t_start_rel_testing + 1
      X <- X_testing_held_out_zeroed[t_start_rel_testing:(S_ - 1), -i, drop = FALSE]
    } else {
      T_ <- intervals_include[k, "T_", drop = T]
      X <- X_testing_held_out_zeroed[t_start_rel_testing:(t_end_rel_testing - 1), -i, drop = FALSE]
    }
    rel_testing_range <- t_start_rel_testing:min(t_end_rel_testing - 1, S_ - 1)
    
    ## comment
    non_z_cum_input_vec_no_held_out_contrib <- c(0, scaled_sum_mat_mult(delta, X, estimated_parameters$W[, i])) +
      estimated_parameters$eta[i] * (1 - delta ^ (1:T_)) / (1 - delta)
    
    #################################################################
    # GET RELATIVE TIMES s FOR LIKELIHOOD COMPS ON CURRENT INTERVAL #
    #################################################################
    
    s_df_cur <- s_intervals_df %>% filter(index == interval_index)
    msg <- paste0(k, " - T = ", T_, " - s rows = ", nrow(s_df_cur))
    print(msg)
    write(paste(msg), update_filename_w_path, append = TRUE)
    
    s_frwrd_incl <- s_df_cur %>% filter(include_frwrd) %>% pull(s)
    s_bckwrd_incl <- s_df_cur %>% filter(include_bckwrd) %>% pull(s)
    # note: whenever s incl in interval, this interval will be backward included
    s_in_interval <- s_df_cur %>% filter(partial_wrt_s) %>% pull(s)
    
    # has to be at least one based on construction of s_intervals_df
    RUN_BCKWRD <- length(s_bckwrd_incl) > 0
    RUN_FRWRD <- length(s_frwrd_incl) > 0
    
    #############################
    # SET IMPUTED HELD OUT COLS #
    #############################
    
    if (RUN_BCKWRD){
      next_row_bckwrds_store <- sum(!is.na(bckwrd_scores_mat_store[,1])) + 1
      bckwrd_scores_mat_store_cur <-
        bckwrd_scores_mat_store[next_row_bckwrds_store + seq_along(s_bckwrd_incl) - 1, , drop = FALSE]
      bckwrd_scores_mat_store_cur[, "index"] <- interval_index
      bckwrd_scores_mat_store_cur[, "s"] <- s_bckwrd_incl
      bckwrd_scores_mat_store_cur[, "partial"] <-
        s_df_cur %>% filter(include_bckwrd) %>% pull(partial_wrt_s) %>% ifelse(1, 0)
      
      B1_mat <- t(sims$COND_PROB_1_MAT[s_bckwrd_incl, rel_testing_range, drop = F])
      colnames(B1_mat) <- paste0("B1_", s_bckwrd_incl)
      B0_mat <- t(sims$COND_PROB_0_MAT[s_bckwrd_incl, rel_testing_range, drop = F])
      colnames(B0_mat) <- paste0("B0_", s_bckwrd_incl)
    } else {
      B1_mat <- NULL
      B0_mat <- NULL
    }
    
    if (RUN_FRWRD){
      if(s_frwrd_incl[1] == 0 & length(s_frwrd_incl) > 1){
        F1_mat <- cbind(0, t(sims$COND_PROB_1_MAT[s_frwrd_incl[s_frwrd_incl > 0], rel_testing_range, drop = F]) *
                          sapply(s_frwrd_incl[s_frwrd_incl > 0], function(s) rel_testing_range <= s))
      } else if (s_frwrd_incl[1] == 0) {
        F1_mat <- matrix(0, nrow = T_ - 1, ncol = 1)
      } else {
        F1_mat <- t(sims$COND_PROB_1_MAT[s_frwrd_incl, rel_testing_range, drop = F]) *
          sapply(s_frwrd_incl, function(s) rel_testing_range <= s)
      }
      colnames(F1_mat) <- paste0("F1_", s_frwrd_incl)
      
      U_col_mat <- matrix(sims$PROB_SPIKE_VEC[rel_testing_range], ncol = 1)
      colnames(U_col_mat) <- "U"
    } else {
      F1_mat <- NULL
      U_col_mat <- NULL
    }
    
    # cols to include in frwrd calcs (also to calc bckwrd partial probs / right censored probs)
    imputed_cols_mat <- cbind(U_col_mat, F1_mat, B1_mat, B0_mat)
      
    # scale cols, add non-held out input
    non_z_cum_input_columns <-
      non_z_cum_input_vec_no_held_out_contrib +
      estimated_parameters$W[get_source_index_in_W(held_out_node, i), i] *
      rbind(
        0, scaled_sum_mat(delta, imputed_cols_mat)
      )
    
    # efficiently draw from non-anticipating distribution
    # to estimate spike probabilities
    no_spike_yet_MAT_sim <- sim_no_spikes_yet_SAPPLY(
      n_sims = n_sims,
      delta = delta,
      kappa = kappa,
      sigma = estimated_parameters$phi["sigma", i],
      non_z_cum_input = non_z_cum_input_columns
    )
    colnames(no_spike_yet_MAT_sim) <- colnames(non_z_cum_input_columns)
    
    
    if (RUN_BCKWRD){
      if (right_censored) {
        #### add comments
        bckwrd_lik_VEC <- no_spike_yet_MAT_sim[T_, c(colnames(B1_mat), colnames(B0_mat))]
      } else {
        #### add comments
        bckwrd_lik_VEC <- no_spike_yet_MAT_sim[T_ - 1,  c(colnames(B1_mat), colnames(B0_mat))] -
          no_spike_yet_MAT_sim[T_,  c(colnames(B1_mat), colnames(B0_mat))]
      }
      
      #### get partial probabilities for intervals that contain s
      ## P(no spike before s) needed for conditional expectation (in denominator)
      ## for interval containing s we scale the stored bckwrd scores (log scale) by this prob
      bckwrd_scores_mat_store_cur[, c("p1", "p0")] <- sim_partial_probs <-
        mapply(
          function(s, partial){
            if (partial == 0){
              return(c(NA, NA))
            } else {
              length_to_check <- s - t_start_rel_testing
              return(no_spike_yet_MAT_sim[length_to_check, c(paste0("B1_", s), paste0("B0_", s))])
            }
          },
          s_bckwrd_incl,  bckwrd_scores_mat_store_cur[, "partial"]
        ) %>% t()

      bckwrd_scores_mat_store_cur[, "lik1"] <- lik1 <- bckwrd_lik_VEC[colnames(B1_mat)]
      bckwrd_scores_mat_store_cur[, "lik0"] <- lik0 <-  bckwrd_lik_VEC[colnames(B0_mat)]
    
      # store raw / initial /old scores (unused but kept for now)
      ## remove later...
      s1 <- log_epsilon_if_small_else_x(lik1, epsilon_0) ->
        bckwrd_scores_mat_store_cur[, "s1"]
      s0 <- log_epsilon_if_small_else_x(lik0, epsilon_0) ->
        bckwrd_scores_mat_store_cur[, "s0"]
      
      ### if probabilities of observing the interval are very small in both cases, discard scores (set both to zero)
      ### make adjustment to very small probabilities
      ### and return score as log of potentially adjusted probabilities
      bckwrd_scores_mat_store_cur[, c("score1", "score0")] <-
        apply(cbind(lik1, lik0, sim_partial_probs), 1,
              function(row) get_adj_score(row, n_sims, sim_ob_threshold, adj_probs)) %>% t()
      
      # store bckwrd scores  
      bckwrd_scores_mat_store[next_row_bckwrds_store + seq_along(s_bckwrd_incl) - 1, ] <-
        bckwrd_scores_mat_store_cur
    }
    
    if (RUN_FRWRD) {
      
      if (right_censored) {
        lik_at_T_VEC <- no_spike_yet_MAT_sim[T_, c("U", colnames(F1_mat))]
      } else {
        lik_at_T_VEC = no_spike_yet_MAT_sim[T_ - 1, c("U", colnames(F1_mat))] -
          no_spike_yet_MAT_sim[T_, c("U", colnames(F1_mat))]
      }
      
      sim_partial_probs_FWRD <-
        mapply(
          function(s, length_to_check){
            if (is.na(length_to_check)){
              return(c(NA, NA))
            } else {
              return(no_spike_yet_MAT_sim[length_to_check, c("U", paste0("F1_", s))])
            }
          },
          s_frwrd_incl,
          vapply(s_frwrd_incl, function(s) ifelse(t_start_rel_testing  < s, s - t_start_rel_testing, NA), 0)
        ) %>% t()
      
      
      for(l in seq_along(s_frwrd_incl)){
        s <- s_frwrd_incl[l]
        #no_spike_yet_s <- no_spike_yet_MAT_sim[, c("U", paste0("F1_", s))]
        no_spike_yet_score_diffs <- no_spike_yet_MAT_sim[, c("U", paste0("F1_", s))] %>%
          apply(1, function(row) get_adj_score(c(row, sim_partial_probs_FWRD[l, ]),
                                               n_sims, sim_ob_threshold, adj_probs,
                                               ret_diff = TRUE))
        
        #lik_at_T_s <- lik_at_T_VEC[c("U", paste0("F1_", s))]
        lik_at_T_score_diff <-
          get_adj_score(c(lik_at_T_VEC[c("U", paste0("F1_", s))], sim_partial_probs_FWRD[l, ]),
                        n_sims, sim_ob_threshold, adj_probs, ret_diff = TRUE)
        
        rs_list_entry <- s_rs_list[[s + 1]]
        rs_incl <- rs_list_entry$rs[t_start_rel_testing < rs_list_entry$rs - 2]
        lengths_include = pmin(t_end_rel_testing, rs_incl - 1) - t_start_rel_testing + 1
        
        add_r_fwrd <- ifelse(
          # check if partial
          t_end_rel_testing >= rs_incl,
          # partial
          no_spike_yet_score_diffs[lengths_include],
          # full
          lik_at_T_score_diff
        )
        
        s_rs_list[[s + 1]]$frwrd_score_sim[paste(rs_incl)] <-
          s_rs_list[[s + 1]]$frwrd_score_sim[paste(rs_incl)] + add_r_fwrd
        
        if (rs_list_entry$p_next_spike_after_cutoff > 0){
          s_rs_list[[s + 1]]$score_next_spike_after_cutoff <-
            s_rs_list[[s + 1]]$score_next_spike_after_cutoff + lik_at_T_score_diff
        }
      }
    }
    
    print(Sys.time() - storeTime)
  }
  print(Sys.time() - startTime)
  
  write("done with intervals - pp computations", update_filename_w_path, append = TRUE)
  print("summing over s")
  s_to_r_frwrd_df <-
    lapply(
      s_vec[!is.na(s_r_max_df$r_max)],
      function(s){
        rs_list_entry <- s_rs_list[[s + 1]]
        s_to_r_probs_unnormalized <-
          exp(rs_list_entry$frwrd_score_sim) * rs_list_entry$s_to_r_NA
          ### CHANGE -- account for times when NA draws end after testing cutoff!
        s_to_r_probs <-
          s_to_r_probs_unnormalized /
          (sum(s_to_r_probs_unnormalized) +
                exp(rs_list_entry$p_next_spike_after_cutoff) * rs_list_entry$p_next_spike_after_cutoff)
        return(
          tibble(
            held_out_node = held_out_node,
            s = s,
            r = rs_list_entry$rs,
            s_to_r_NA = rs_list_entry$s_to_r_NA,
            s_to_r_prob = unname(s_to_r_probs)
          )
        )
      }
    ) %>%
    bind_rows()
  
  print("calculate forward probabilities")
  pp_fwrd_0 <- numeric(length(testing_range) + 1)
  pp_fwrd_0[1] <- 1
  for(r_cur in 1:length(testing_range)){
    s_matches <-
      s_to_r_frwrd_df %>% filter(r == r_cur)
    if (nrow(s_matches) == 0){
      pp_fwrd_0[r_cur + 1] <- 0
    } else {
      pp_fwrd_0[r_cur + 1] <- sum(s_matches$s_to_r_prob * pp_fwrd_0[s_matches$s + 1])
    }
  }
  pp_fwrd <- pp_fwrd_0[-1]
  
  print("calculate estimated (frwrd-bckwrd) probabilities")
  s_intervals_bckwrd_sum <-
    s_intervals_df %>%
    filter(include_bckwrd) %>%
    left_join(
      as_tibble(bckwrd_scores_mat_store),
      by = c("index", "s")
    ) %>%
    group_by(s) %>%
    summarize(
      across(starts_with("score"), sum, .names = "{.col}")
    )
  
  # calculate smoothed posterior predictive probabilities
  pp_vec <- numeric(length(testing_range))
  for (k in seq_along(s_r_max_df$s)[-1]){
    t <- s_r_max_df$s[k]
    pp_fwrd_cur <- pp_fwrd[t]
    
    if(s_r_max_df$p_NA[k] == 0){
      pp_vec[t] <- 0
    } else if (!(t %in% s_intervals_bckwrd_sum$s)) {
      pp_vec[t] <- pp_fwrd_cur
    } else {
      bckwrd_scores <- s_intervals_bckwrd_sum %>% filter(s == t)
      N0 <- bckwrd_scores[1, "score0", drop = T]
      N1 <- bckwrd_scores[1, "score1", drop = T]
      log_pp <- N1 + log(pp_fwrd_cur) - log_x_plus_y(N0 + log(1 - pp_fwrd_cur), N1 + log(pp_fwrd_cur))
      pp_vec[t] <- unname(exp(log_pp))
    }
  }
  
  df_out <- tibble(
    s_r_max_df %>% filter(s > 0) %>% rename(r = s)
  ) %>%
    mutate(
      t = testing_range,
      held_out_node = held_out_node,
      look_ahead_limit_bckwrd = look_ahead_limit_bckwrd,
      n_sims = n_sims,
      frwrd_sim = pp_fwrd,
      pp_sim = pp_vec,
    ) 
  
  NA_sum <- sum(s_r_max_df$p_NA)
  pp_sum <- sum(pp_vec)
  if (pp_sum > NA_sum) {
    df_out$pp_scaled_sim <- pp_vec * NA_sum / pp_sum
  } else {
    df_out$pp_scaled_sim <- pp_vec * (S_ - NA_sum) / (S_ - pp_sum)
  }
  
  if (!is.null(T_train_vec)) {
    df_out <- df_out %>%
      mutate(
        bs_prob = sapply(
          1:10000,
          function(n_sims_full){
            T_s_draw <- sample(x = T_train_vec, size = S_, replace = T)
            spike_times <- cumsum(T_s_draw)[cumsum(T_s_draw) <= S_]
            draw <- numeric(S_)
            draw[spike_times] <- 1
            return(draw)
          }
        ) %>% rowMeans()
      )
  }
  
  if(save_to_file){
    if (is.null(filename_w_path)) {
      saveRDS(
        df_out, paste0("pp_", held_out_node, "_",
                       str_extract_all(Sys.time(), "[0-9]") %>% unlist() %>% paste(collapse = '') %>% substr(1, 8),
                       ".RData")
      )
    } else {
      saveRDS(df_out, filename_w_path)
    }
  }
  
  write(paste("ended at", Sys.time()), update_filename_w_path, append = TRUE)
  
  return(
    list(
      run_time = Sys.time() - startTime,
      df_out = df_out
    )
  )
}

run_frwrd_bckwrd_nodes <- function(
    nodes_to_run = NULL,
    NA_sims_list,
    T_train_cutoff,
    T_test_cutoff,
    #s_vec_mats_list,
    delta = .975,
    kappa = 50,
    n_sims = NULL,
    look_ahead_limit_bckwrd = Inf,
    epsilon_0 = 1e-06, # to prevent zero/neg infinity likelihood calculations, should be no bigger than 1 / nsim...
    epsilon_1 = 1e-25, # for final estimates to enable proper entropy loss calculations
    parallelize_over_nodes = TRUE,
    n_nodes_par = NULL,
    save_to_file = TRUE,
    pp_out_dir = getwd(),
    file_suffix = str_extract_all(Sys.time(), "[0-9]") %>% unlist() %>% paste(collapse = '') %>% substr(1, 8),
    sim_ob_threshold = 10,
    upper_tail_prob_for_adj = .65,
    T_LB = 0,
    T_UB = Inf
  ){
  
  if (!("updateProgress" %in% list.files(pp_out_dir))){
    dir.create(paste0(pp_out_dir, "/updateProgress"))
  }
  
  if(!("filter_length" %in% names(NA_sims_list$spike_intervals_df))){
    NA_sims_list$spike_intervals_df$filter_length <- "keep"
  }
  
  estimated_parameters <- NA_sims_list$estimated_parameters
  
  #last_observed_spike_times <- list_with_held_out_sims$last_observed_spike_times
    #setNames(list_with_held_out_sims$last_observed_spike_times$last_spike,
    #         list_with_held_out_sims$last_observed_spike_times$i)
  
  
  T_train_list <- mapply(
    function(node_i, last_observed_spike_time){
      NA_sims_list$spike_intervals_df %>%
        filter(t_end < last_observed_spike_time & i == node_i) %>%
        pull(T_)
    },
    nodes_to_run,
    vapply(nodes_to_run, function(node_i) NA_sims_list[[paste(node_i)]]$last_observed_spike_time, 0),
    SIMPLIFY = F
  )
  names(T_train_list) <- paste(nodes_to_run)
  
  if(parallelize_over_nodes){
    cl <- makeCluster(n_nodes_par, type="FORK")
  } else {
    cl <- NULL
  }
  outList <-
    lapply_parLapplyLB(
      cl,
      seq_along(nodes_to_run),
      function(k){
        held_out_node <- nodes_to_run[k]
        node_str <- paste(held_out_node)
        testing_range <- NA_sims_list[[node_str]]$testing_range
        last_observed_spike_time <- NA_sims_list[[node_str]]$last_observed_spike_time
        
        get_fwrd_bckwrd(
          held_out_node = held_out_node,
          sims = NA_sims_list[[node_str]]$NA_sims,
          testing_range = testing_range,
          estimated_parameters = estimated_parameters,
          #s_vec_mats_list = s_vec_mats_list,
          intervals_testing_excl_held_out =
            NA_sims_list$spike_intervals_df %>%
            filter(
              t_start > last_observed_spike_time &
                t_start < T_test_cutoff - 1 &
                i != held_out_node &
                filter_length == "keep" &
                T_ > T_LB & T_ < T_UB
            ) %>%
            mutate(
              t_start_testing_period = t_start - last_observed_spike_time,
              t_end_testing_period = t_end - last_observed_spike_time,
              interval_partially_observed = t_end > T_test_cutoff
            ),
          binary_spike_mat_testing_period = NA_sims_list$full_spike_mat[testing_range, ],
          delta = delta,
          kappa = kappa,
          n_sims = n_sims,
          look_ahead_limit_bckwrd = look_ahead_limit_bckwrd,
          epsilon_0 = epsilon_0,
          epsilon_1 = epsilon_1,
          T_train_vec = T_train_list[[paste(held_out_node)]],
          save_to_file = save_to_file,
          filename_w_path = paste0(pp_out_dir, "pp_", held_out_node, "_", file_suffix, ".RData"),
          update_filename_w_path = paste0(pp_out_dir, "/updateProgress/", node_str, ".txt"),
          sim_ob_threshold = sim_ob_threshold,
          upper_tail_prob_for_adj = upper_tail_prob_for_adj
        )
      }
    )
  if(parallelize_over_nodes){
    stopCluster(cl)
  }
  
  return(outList)
}

# function for taking the logarithm of 2 very small numbers,
# which are passed to the function logged
log_x_plus_y <- function(log_x, log_y){
  max(log_x, log_y) + log1p(exp(-abs(log_x - log_y)))
}
