log_joint_cond <- function(z, delta, non_z_cum_input_vec, sigma, kappa){
  T_ <- length(z)
  z_scaled_sums <- vapply(1:T_, function(ell) scaled_sum(rate = delta, vec = z[1:ell]), 0)
  
  v <- non_z_cum_input_vec + z_scaled_sums
  
  ll <- -T_ * log(sigma) - .5 * (1 / sigma ^ 2) * sum(z ^ 2) -
    sum(log(1 + exp(kappa * (v - 1)))) + kappa * v[T_]
  return(ll)
}
log_join_lik_interval <- function(z, v, T_, sigma, kappa = 50, delta = .975){
  dnorm(z, mean = 0, sd = sigma, log = T) - sum(log(1 + exp(kappa * (v - 1)))) + kappa * v[T_]
}
rejection_sample_try <- function(T_, delta, non_z_cum_input_vec, sigma, kappa){
  
  #candidate
  z <- rnorm(T_, mean = 0, sd = sigma)
  
  # get v
  z_scaled_sums <- vapply(1:T_, function(ell) scaled_sum(rate = delta, vec = z[1:ell]), 0)
  v <- non_z_cum_input_vec + z_scaled_sums
  
  # uniform sample
  u <- runif(1)
  
  accept <-
    log(u) <= kappa * v[T_] - sum(log(1 + exp(kappa * (v - 1)))) - kappa
  
  return(
    list(
      z=z,
      v=v,
      ll=log_joint_cond(z, delta, non_z_cum_input_vec, sigma, kappa),
      accept=accept
    )
  )
}


#sim = T
#interval_index <- 102
#interval_index <- 29


eval_interval_sampling <- function(
  interval_index,
  runList,
  n_sims_q = 1000,
  get_AR_samples = F,
  n_proposals_AR = 1000,
  n_sims_AR = 1000,
  num_tries = 200,
  sim = TRUE,
  return_samples = F,
  draw_from_q0 = F,
  #empirical_KL = F,
  AR_lower_lim = 0,
  get_truth = T,
  only_for_PSIS_calc = F
){
  
  return_list <- list()
  
  if (sim) {
    # get observable data from given interspike interval
    interval_data <- get_interval_data(interval_index, spike_intervals = runList$model_params$sim$spike_intervals,
                                       binary_spike_train = runList$model_params$sim$binary_spike_train)
    return_list$i <- i <- interval_data$i
    return_list$T_ <- T_ <- interval_data$T_
    X <- interval_data$X
    
    if (get_truth){
      #################
      # get actual latent path data for this interval
      start.end <- runList$model_params$sim$spike_intervals[interval_index, c("t_start", "t_end")] %>% as.numeric()
      actual <- runList$model_params$sim$outMatList[[i]][start.end[1]:start.end[2], ]
      V_T <- actual[T_, "V"]
      return_list$z_true <- z_true <- actual[, "z"]
      return_list$v_true <- v_true <- actual[, "V"]
      # cumumlative scaled z sums
      #z_scaled_sums <- vapply(1:T_, function(ell) scaled_sum(rate = delta, vec = z[1:ell]), 0)
      #ll_actual <- log_joint_cond(actual[, "z"], delta, non_z_cum_input_vec, sigma, kappa)
      ################## 
    }
    
    ### KNOWN MODEL PARAMS
    delta <- runList$model_params$delta
    kappa <- runList$model_params$K
    
  } else {
    # get observable data from given interspike interval
    interval_data <- get_interval_data(interval_index, spike_intervals = runList$fit_params$spike_intervals_df_train,
                                       binary_spike_train = runList$fit_params$binary_spike_mat_train)
    return_list$i <- i <- interval_data$i
    return_list$T_ <- T_ <- interval_data$T_
    X <- interval_data$X
    
    ### MODEL PARAMS ASSUMED KNOWN
    delta <- runList$fit_params$delta
    kappa <- runList$fit_params$K
  }
  
  ### ESTIMATED MODEL PARAMETERS
  W_in <- runList$adam_update$values_out$W[, i]
  eta <- runList$adam_update$values_out$eta[i]
  phi <- runList$adam_update$values_out$phi[, i]
  sigma <- phi["sigma"]
  
  q_fit_list <- obj_grad_single_spike(W_in, phi, eta, T_, X, delta = .975, K=50, grads = F, retFittedVals = T)
  
  return_list$q_fit_list <- q_fit_list
  
  non_z_cum_input_vec <- c(0, scaled_sum_mat_mult(delta, X, W_in)) +
    eta * (1 - delta ^ (1:T_)) / (1 - delta)
  
  # simulate from q dist
  Z_sim_q <- mat.mult(q_fit_list$R_out$R, matrnorm(T_, n_sims_q)) + q_fit_list$mus
  scaled_sum_z_sim_q <- scaled_sum_mat(delta, Z_sim_q)
  V_sim_q <- non_z_cum_input_vec + scaled_sum_z_sim_q
  
  return_list$log_lik_mat <- cbind(
    # log likelihood of Z_sim_q under q
    "log_q_Z_sim_q" =
      mvnfast::dmvn(t(Z_sim_q), q_fit_list$mus, q_fit_list$R_out$R, log = T, ncores = 1, isChol = TRUE),
    # log joint likelihood of Z_sim_q and y
    "log_joint_Z_sim_z" =
      colsums(dnorm(Z_sim_q, mean = 0, sd = sigma, log = T)) -
      colsums(Log(1 + exp(kappa * (V_sim_q - 1)))) +
      kappa * V_sim_q[T_, ]
  )
  
  if (only_for_PSIS_calc){
    return(return_list$log_lik_mat)
  }
  
  if(return_samples){
    #matices
    return_list$Z_sim_q <- Z_sim_q
    return_list$V_sim_q <- V_sim_q
  }
  
  if (sim) {
    # log likelihood of true Z under q
    return_list$log_q_Z_actual <-
      mvnfast::dmvn(return_list$z_true, q_fit_list$mus, q_fit_list$R_out$R, log = T, ncores = 1, isChol = TRUE)
    # log joint likelihood of true Z and y
    return_list$log_joint_Z_sim_actual <-
      sum(dnorm(return_list$z_true, mean = 0, sd = sigma, log = T)) -
      sum(log(1 + exp(kappa * (return_list$v_true - 1)))) +
      kappa * return_list$v_true[T_]
    return_list$mahala_true_qz <- 
      mahala(matrix(return_list$z_true, nrow = 1), q_fit_list$mus, sigma = t(q_fit_list$R_out$R), ischol = T)
    return_list$mahala_true_Pz <- 
      mahala(matrix(return_list$z_true, nrow = 1), numeric(T_), sigma = sigma ^ 2 * diag(T_), ischol = F)
  }
  
  if (get_AR_samples){
    #ar_return <- list()
    accepted_zs_mat <- matrix(NA, nrow = T_, ncol = n_sims_AR)
    accepted_vs_mat <- matrix(NA, nrow = T_, ncol = n_sims_AR)
    sims_remaining <- n_sims_AR
    #success <- FALSE
    
    return_list$abort_AR = F
    for (make_try in 1:num_tries){
      z_mat_try <-  matrix(Rnorm(n_proposals_AR * T_, m = 0, s = sigma), nrow = T_)
      scaled_sum_z_mat_try <- scaled_sum_mat(delta, z_mat_try)
      v_mat_try <- non_z_cum_input_vec + scaled_sum_z_mat_try
      
      # uniform sample
      u <- runif(n_proposals_AR)
      which_accept <- which(
        log(u) + kappa <= kappa * v_mat_try[T_, ] - colsums(Log(1 + exp(kappa * (v_mat_try - 1))))
      )
      num_accept <- length(which_accept)
      #cat(num_accept, "accepted out of", n_sims_try, "draws,", min(sims_remaining - num_accept), "remaining\n")
      if (num_accept > 0 & num_accept > AR_lower_lim) {
        if (num_accept >= sims_remaining){
          accepted_zs_mat[, (n_sims_AR - sims_remaining + 1):n_sims_AR] <-
            columns(z_mat_try, which_accept)[, 1:sims_remaining]
          
          accepted_vs_mat[, (n_sims_AR - sims_remaining + 1):n_sims_AR] <-
            columns(v_mat_try, which_accept)[, 1:sims_remaining]
          
          sims_remaining <- 0
          #success <- TRUE
          break
        } else {
          accepted_zs_mat[, n_sims_AR - sims_remaining + 1:num_accept] <-
            columns(z_mat_try, which_accept)
          accepted_vs_mat[, n_sims_AR - sims_remaining + 1:num_accept] <-
            columns(v_mat_try, which_accept)
          sims_remaining <- sims_remaining - num_accept
        }
      } else {
        return_list$abort_AR = T
        print(paste("less than", AR_lower_lim + 1, "accepts in", n_proposals_AR, "draws, ABORT!"))
        break
      }
    }
    return_list$AR_sims_returned <- AR_sims_returned <- n_sims_AR - sims_remaining
    return_list$tries <- make_try
    
    if (AR_sims_returned > 0) {
      Z_sim_AR <- columns(accepted_zs_mat, 1:AR_sims_returned)
      V_sim_AR <- columns(accepted_vs_mat, 1:AR_sims_returned)
      
      if(return_samples){
        #matices
        return_list$Z_sim_AR <- Z_sim_AR
        return_list$V_sim_AR <- V_sim_AR
      }
      
      # mahalanobis distances to q, comparisons
      return_list$maha_mat <- cbind(
        "samp_ar_wrt_q" =
          mahala(transpose(Z_sim_AR), q_fit_list$mus, sigma = t(q_fit_list$R_out$R), ischol = T),
        "samp_q_wrt_q" =
          mahala(transpose(columns(Z_sim_q, 1:AR_sims_returned)), q_fit_list$mus, sigma = t(q_fit_list$R_out$R), ischol = T),
        "samp_ar_wrt_pZ" =
          mahala(transpose(Z_sim_AR), numeric(T_), sigma = sigma ^ 2 * diag(T_), ischol = F)
      )
      
      return_list$v_star_RA <- .5 * (V_sim_AR[T_,] + V_sim_AR[T_ - 1])
      
      return_list$log_lik_mat_AR <- cbind(
        # log likelihood of Z_sim_AR under q
        "log_q_Z_sim_AR" =
          mvnfast::dmvn(t(Z_sim_AR), q_fit_list$mus, q_fit_list$R_out$R, log = T, ncores = 1, isChol = TRUE),
        # log joint likelihood of Z_sim_AR and y
        "log_joint_Z_sim_AR" =
          colsums(dnorm(Z_sim_AR, mean = 0, sd = sigma, log = T)) -
          colsums(Log(1 + exp(kappa * (V_sim_AR - 1)))) +
          kappa * V_sim_AR[T_, ]
      )
      
      #if (empirical_KL) {
      #  Z_sim_0 <- matrix(Rnorm(n = T_ * n_sims_q, m = 0, s = sigma), nrow = T_)
      #  return_list$emirical_KL <-
      #    cbind(
      #      "samp_q_rel_samp_AR" = FNN::KL.divergence(t(Z_sim_q), t(Z_sim_AR)),
      #      "samp_pZ_rel_samp_AR" = FNN::KL.divergence(t(Z_sim_0), t(Z_sim_AR))
      #    )
      #}
    }
  }
  if(draw_from_q0) {
    # simulate from vanilla z ~ N(0, sigma) for comparison
    return_list$Z_sim_0 <- matrix(Rnorm(n = T_ * n_sims_q, m = 0, s = sigma), nrow = T_)
    return_list$V_sim_0 <- scaled_sum_mat(delta, return_list$Z_sim_0) + non_z_cum_input_vec
  }
  
  return(return_list)
}


#cur_interval_index
#df_w_psis_out_ROW <- training_intervals_df_w_psis_out[cur_interval_index,]

get_plots <- function(cur_interval_index, df_w_psis_out_ROW, chi_sq_label_x_coord = 0, scale = 1, snr_rev_rank = 1:100){
  
  interval_out <- 
    eval_interval_sampling(
      interval_index = cur_interval_index,
      runList = runList_50000,
      n_sims_q = 10000,
      get_AR_samples = T,
      n_proposals_AR = 5000,
      n_sims_AR = 10000,
      num_tries = 1000,
      sim = TRUE,
      return_samples = T,
      draw_from_q0 = T,
      #empirical_KL = F,
      AR_lower_lim = 1
    )
  
  pMaha <- tibble(
    mahalanobis = interval_out$maha_mat[, "samp_ar_wrt_q"]
  ) %>%
    ggplot(mapping = aes(x = mahalanobis)) +
    geom_histogram(
      #breaks = seq(from = plot_x_lb, to = plot_x_ub, by = .02),
      mapping = aes(y = after_stat(density)),
      fill="cornflowerblue"
    ) +
    stat_function(fun = function(x) dchisq(x, df = interval_out$T_),
                  color = "darkgoldenrod1", linewidth = .75, alpha = .9, n = 301) +
    geom_vline(
      xintercept = interval_out$mahala_true_qz,
      color = "firebrick1"
    ) +
    scale_y_continuous(expand = c(.01,0)) +
    scale_x_continuous(expand = c(0,0),
                       breaks = c(0, 25, 50, 75, 100, 150, 175, 200),
                       labels = c("0", "", "", "75", "", "", "150", "")) +
    labs(
      #subtitle = TeX(r'(Distances between $q_{\hat{\phi}}(z)$ and draw from $p_{\hat{\theta},\hat{\bf{W}}}(z|y)$)'),
      subtitle = TeX(paste0(#"Inter-spike period (i=",
                         "(i=",
                         snr_rev_rank[df_w_psis_out_ROW[1, "i", drop = T]],", k=", df_w_psis_out_ROW[1, "k", drop = T],
                         "), $\\Delta t_k^i$ = ", df_w_psis_out_ROW[1, "T_", drop = T], ", $\\hat{k}$ = ",
                         round(df_w_psis_out_ROW[1, "pareto_k_interval", drop = T],2))),
      #title = "Mahalanobis distances",
      y = "",
      x = ""
      ) +
    theme_bw() +
    theme(
      #plot.title = element_text(hjust = .5),
      plot.subtitle = element_text(hjust = .5),
      plot.title = element_blank(),
      plot.margin = margin(0, 0, 0, 0, "cm")
    ) +
    geom_text(
      data = tibble(label = "distance to\ntrue z", x = interval_out$mahala_true_qz - 0 * interval_out$T_ / 9, y = .01),
      mapping = aes(x = x, y = y, label = label), color = "black"#, hjust = "inward", vjust = "inward"
    ) +
    annotate(
      "text", x=chi_sq_label_x_coord, y=.001,
      label = TeX(paste0("$\\chi_{", interval_out$T_, "}$")),
      color = "darkgoldenrod1", hjust = "inward", vjust = "inward",
      size = 5
    )
  
  sim_labels_z <- c('z~"~"~p[list(hat(theta),hat(W))]*(z~"|"~Y)',
                    'z~"~"~q[hat(phi)]*(z)')

  sim_labels_v <- c('z~"~"~p[list(hat(theta),hat(W))]*(z~"|"~Y) %->% v',
                  'z~"~"~q[hat(phi)]*(z) %->% v')
  
  mean_q_labels_v <- c("mean~AR~draws", 'E[q[hat(phi)]]*(v)')
  mean_q_labels_z <- c("mean~AR~draws", 'E[q[hat(phi)]]*(z)')
  
  df_w_sims <- 
    setNames(as_tibble(interval_out$V_sim_AR[,1:100]), paste0("draw", 1:100)) %>%
    mutate(
      source = "AR_samp_v",
      latent_var = "v",
      t=row_number()
    ) %>%
    bind_rows(
      setNames(as_tibble(interval_out$V_sim_q[,1:100]), paste0("draw", 1:100)) %>%
        mutate(
          source = "sim_q_v",
          latent_var = "v",
          t=row_number()
        ),
      setNames(as_tibble(interval_out$Z_sim_AR[,1:100]), paste0("draw", 1:100)) %>%
        mutate(
          source = "AR_samp_z",
          latent_var = "z",
          t=row_number()
        ),
      setNames(as_tibble(interval_out$Z_sim_q[,1:100]), paste0("draw", 1:100)) %>%
        mutate(
          source = "sim_q_z",
          latent_var = "z",
          t=row_number()
        )
    ) %>%
    pivot_longer(cols = 1:100, names_to = "draw", values_to = "sampled_latent_var") %>%
    mutate(
      source = factor(source, levels = c("AR_samp_z", "sim_q_z", "AR_samp_v", "sim_q_v"),
                      labels = c(sim_labels_z, sim_labels_v))
      #labels = c('z~"~"~p[list(hat(theta),hat(W))]*(z~"|"~Y) %->% v',
      #           'z~"~"~q[hat(phi)]*(z) %->% v'))
    )
  
  p_z <- 
    df_w_sims %>%
    filter(latent_var == "z") %>%
    ggplot(mapping = aes(x = t, y = sampled_latent_var, group = draw)) +
    facet_wrap(~source, labeller = label_parsed) +
    geom_line(alpha = .4) +
    geom_line(
      inherit.aes = F,
      data =
        tibble(sim_label1 = sim_labels_z[1],
               sim_label2 = sim_labels_z[2],
               t = 1:interval_out$T_,
               sampled_latent_var=interval_out$z_true) %>%
        pivot_longer(cols = 1:2, names_to = "names",
                     values_to = "source"),
      mapping = aes(x=t, y=sampled_latent_var),
      color = "firebrick1",
      linewidth = .5
    ) +
    geom_line(
      inherit.aes = F,
      data =
        bind_rows(
          tibble(source = sim_labels_z[1],
                 t = 1:interval_out$T_,
                 mean_z = rowMeans(interval_out$Z_sim_AR)),
          tibble(source = sim_labels_z[2],
                 t = 1:interval_out$T_,
                 mean_z = interval_out$q_fit_list$mus)
        ),
      mapping = aes(x=t, y=mean_z),
      color = "dodgerblue",
      linewidth = .5
    ) +
    geom_text(
      inherit.aes = F,
      data = tibble(source = sim_labels_z,
                    label = mean_q_labels_z,
                    x = interval_out$T_ - .5,
                    y = min(filter(df_w_sims, latent_var == "z")$sampled_latent_var, interval_out$z_true) - .02 * scale),
      parse = T,
      mapping = aes(x=x, y=y, label = label),
      color = "dodgerblue",
      vjust = "inward",
      hjust = "inward"
    ) +
    geom_text(
      inherit.aes = F,
      data = tibble(source = sim_labels_z,
                    label = "true sim data",
                    x = (interval_out$T_ + 1) / 2,
                    y = max(filter(df_w_sims, latent_var == "z")$sampled_latent_var, interval_out$z_true) + .02 * scale),
      mapping = aes(x=x, y=y, label = label),
      color = "firebrick1",
      vjust = "inward",
      hjust = "inward"
    ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_text(hjust = .5),
      plot.margin = margin(.1, 0, 0, .1, "cm")
    ) +
    labs(
      #subtitle = "Sampled latent variables",
      #subtitle = "with true values from the simulated data",
      y = TeX(r'((b) $z_{(i,k)}$ draws)'),
      x = ""
    ) +
    scale_x_continuous(expand=c(0,0))
  
  if(scale == 1){
    p_z <- p_z  + scale_y_continuous(breaks = seq(-.15, .20, by =.05),
                                     labels =c("-.15", "-.10", "-.05", "0.00", "0.05","0.10","0.15", "0.20"),
                                     expand=c(.005 * scale, .005 * scale))
  } else if (scale < 1) {
    p_z <- p_z  + scale_y_continuous(breaks = seq(-.06, .06, by =.02),
                                     labels =c("-.06", "-.04", "-.02", "0.00", "0.02","0.4","0.6"),
                                     expand=c(.005 * scale, .005 * scale))
  } else {
    p_z <- p_z  + scale_y_continuous(breaks = seq(-.20, .20, by =.1),
                                     labels =c("-.20", "-.10", "0.00", "0.10", "0.20"),
                                     expand=c(.005 * scale, .005 * scale))
  }
  
  p_v <- df_w_sims %>%
    filter(latent_var == "v") %>%
    ggplot(mapping = aes(x = t, y = sampled_latent_var, group = draw)) +
    facet_wrap(~source, labeller = label_parsed) +
    geom_line(alpha = .4) +
    geom_line(
      inherit.aes = F,
      data =
        tibble(sim_label1 = sim_labels_v[1],
               sim_label2 = sim_labels_v[2],
               t = 1:interval_out$T_,
               sampled_latent_var=interval_out$v_true) %>%
        pivot_longer(cols = 1:2, names_to = "names",
                     values_to = "source"),
      mapping = aes(x=t, y=sampled_latent_var),
      color = "firebrick1",
      linewidth = .5
    ) +
    geom_line(
      inherit.aes = F,
      data =
        bind_rows(
          tibble(source = sim_labels_v[1],
                 t = 1:interval_out$T_,
                 mean_v = rowMeans(interval_out$V_sim_AR)),
          tibble(source = sim_labels_v[2],
                 t = 1:interval_out$T_,
                 mean_v = interval_out$q_fit_list$m_vec)
        ),
      mapping = aes(x=t, y=mean_v),
      color = "dodgerblue",
      linewidth = .5
    ) +
    geom_text(
      inherit.aes = F,
      data = tibble(source = sim_labels_v,
                 label = mean_q_labels_v,
                 x = interval_out$T_ - .5,
                 y = min(filter(df_w_sims, latent_var == "v")$sampled_latent_var) - .08), #+.03
      parse = T,
      mapping = aes(x=x, y=y, label = label),
      color = "dodgerblue",
      vjust = "inward",
      hjust = "inward"
    ) +
    geom_text(
      inherit.aes = F,
      data = tibble(source = sim_labels_v,
                   label = "true sim data",
                   x = (interval_out$T_ + 1) / 2,
                   y = max(filter(df_w_sims, latent_var == "v")$sampled_latent_var) + .08), #-.01
      mapping = aes(x=x, y=y, label = label),
      color = "firebrick1",
      vjust = "inward",
      hjust = "inward"
    ) +
    theme_bw() +
    theme(
      #plot.title = element_text(hjust = .5),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(0, 0, 0, .1, "cm"),
      #plot.title.position = "none"
    ) +
    labs(
      title = "",
      subtitle = "",
      y = TeX(r'((c) $v_{(i,k)}$ draws)'),
      #x = "t (within interspike interval)",
      x = ""
    ) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(breaks = seq(-0, 1, by =.25), expand=c(.03,.03))

  return(list(pMaha, p_z, p_v))
}


