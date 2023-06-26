root_dir <- "/Users/williamreedpalmer/Documents/Columbia/Research/inferring network VI/R code"
plot_dir <- paste0(root_dir, "/0_output/sim/plots/")
setwd(root_dir)

eta_vec <- rep(c(.025, .03, .035, .04), 5)
sigma_vec <- c(rep(.01, 4), rep(.02, 4), rep(.03, 4), rep(.04, 4), rep(.05, 4))

Ts <- c(1000, 5000, 10000, 25000, 50000)

# estimates get from fit
###
# for each directed pair $i \ne j$ 
# take mean and sd of estimates from final 25 epochs
# determine that the model estimates a postive edges from i to j if
# the mean - 3*sd > 0 and mean - 3*sd > .05
# the second condition accounts for the convergence of estimates very close to zero
# analagously determine that the model estimates a negative edges from i to j if
# the mean + 3*sd < 0 and mean + 3*sd < .05

mean_sd_last_25 <- W_paths %>%
  filter(epochs_run - step < 25) %>%
  group_by(run_name, source_node, target_node, Wij_true) %>%
  summarize(mean_est = mean(Wij_est), sd_est = sd(Wij_est))

# really bad merge
compare_estimated_signs <- mean_sd_last_25 %>%
  mutate(
    ci_lb = mean_est - 3 * sd_est,
    ci_ub = mean_est + 3 * sd_est,
    est_sign = case_when(
      ci_ub < 0 & ci_lb < -.05 ~ -1,
      ci_lb > 0 & ci_ub > .05 ~ 1,
      TRUE ~ 0
    ),
    check = est_sign == sign(Wij_true)
  )
compare_estimated_signs <-
  compare_estimated_signs %>%
  left_join(
    W_paths %>%
      filter(epochs_run - step == -1) %>%
      select(sqrd_error, source_node, target_node, run_name)
  )
VB_out <- compare_estimated_signs


table(compare_estimated_signs$run_name, compare_estimated_signs$check)

#### COMPUTE TEs
# slow so compute only 1k - 25k for now...
get_TE <- function(i,j){
  lapply(
    Ts,
    function(T_train){
      RTransferEntropy::transfer_entropy(
        sim_on_network$binary_spike_train[1:T_train,i],
        sim_on_network$binary_spike_train[1:T_train,j],
        nboot = 100, lx = 1, ly = 1)$coef %>%
        as_tibble() %>%
        mutate(
          T_train = T_train,
          source_node = c(i, j),
          target_node = c(j, i)
        )
    }
  ) %>%
    bind_rows()
}

#### COMPUTE TEs --- SLOOOOOOWW
if (FALSE){
  library(RTransferEntropy)
  i_j_combos <- compare_estimated_signs %>%
    ungroup() %>%
    filter(run_name == "25000") %>%
    select(source_node, target_node)
  
  cl <- makeCluster(13, type="FORK")
  TE_out_ALL <- 
    parLapply(
      cl,
      seq_along(i_j_combos$source_node),
      function(k){
        get_TE(i_j_combos$source_node[k], i_j_combos$target_node[k])
      }) %>%
    bind_rows()
  stopCluster(cl)
  
  saveRDS(TE_out_ALL, "TE_out.RData")
}
TE_out_ALL <- readRDS("TE_out.RData")


## get GLM fits for comparison method
## WARNING -- Long runtime
if (FALSE){
  out_df <- 
    lapply(seq_along(sim_on_network$outMatList),
           function(k){
             sim_on_network$outMatList[[k]] %>% as_tibble %>% mutate(node_i = k)
           }) %>%
    bind_rows()
  
  fit_lists <- list()
  for (i in 1:20){
    print(i)
    fit_list_node <- list()
    for (T_fit in Ts){
      print(T_fit)
      fit_list_node[[paste0("T_", T_fit)]] <- list()
      fit_list_node[[paste0("T_", T_fit)]]$full <-
        glm(
          y ~ .,
          family = binomial(link = "logit"),
          data = out_df %>% filter(node_i == i, t <= T_fit) %>%
            select(y, intercept) %>%
            bind_cols(
              setNames(
                as_tibble(sim_on_network$networkPredictors[[i]][1:T_fit, -i]),
                paste0("n",(1:n)[-i])
              )
            )
        )
      fit_list_node[[paste0("T_", T_fit)]]$stepAIC_out =
        MASS::stepAIC(fit_list_node[[paste0("T_", T_fit)]]$full, trace = F)
    }
    fit_lists[[i]] <- fit_list_node
  }
  saveRDS(fit_lists, "glm_fits.RData")
}
fit_lists <- readRDS("glm_fits.RData")

coef_list <- list()
iter = 0
for (node_i in 1:20){
  for (T_fit in Ts){
    
    for(fit in c("stepAIC_out", "full")){
      iter = iter + 1
      c_vec <- setNames(numeric(n), paste0("n",(1:n)))
      coef <- fit_lists[[node_i]][[paste0("T_", T_fit)]][[fit]]$coefficients
      for(k in seq_along(coef[-c(1,2)])){
        cur_var <- names(coef)[k + 2]
        c_vec[cur_var] <- coef[cur_var]
      }
      coef_list[[iter]] <- as_tibble(t(c_vec)) %>% mutate(node_i = node_i, T_fit = T_fit, model = fit)
    }
  }
}

GLM_results <- bind_rows(coef_list) %>%
  pivot_longer(cols = 1:20, names_to = "source", values_to = "coef") %>%
  mutate(
    coef_scale = coef / 50,
    source_node = as.numeric(str_sub(source, start = 2))
  ) %>%
  rowwise() %>%
  mutate(
    true_val = W_true[source_node, node_i],
    is_edge = !(true_val == 0),
    sqrd_error = (true_val - coef_scale) ^ 2,
    true_sign = sign(true_val),
    correct_sign = sign(coef) == true_sign,
    correct = sign(coef) ^ 2 == true_sign ^ 2,
    target_sigma = sigma_vec[node_i]
  ) %>%
  filter(node_i != source_node)

GLM_results %>%
  group_by(T_fit, model) %>%
  summarize(
    sqrd_error = sum(sqrd_error),
    error_rate_sign = mean(!correct_sign),
    error_rate_incl = mean(!correct)
  )

VB_out %>%
  ungroup() %>%
  filter(run_name == 50000) %>%
  rename(node_i = target_node) %>%
  group_by(node_i) %>%
  summarize(sqrd_error = sum(sqrd_error))

GLM_VB_df_w_SNR <-
  GLM_results %>%
  ungroup() %>%
  filter(T_fit == 50000 & model == "stepAIC_out") %>%
  group_by(node_i, target_sigma) %>%
  summarize(GLM_fit = sum(sqrd_error)) %>%
  ungroup() %>%
  mutate(
    snr_rev_rank = snr_rev_rank_permutation,
    snr_dB = est_out$SNR_db,
    snr_dB_from_fit = est_out$SNR_db_from_fit
  ) %>%
  left_join(
    VB_out %>%
      ungroup() %>%
      filter(run_name == 50000) %>%
      rename(node_i = target_node) %>%
      group_by(node_i) %>%
      summarize(VB_fit = sum(sqrd_error)),
    by = "node_i"
  )

legend_labels = TeX(c("$SNR_i^{\\bf{W},\\bf{Y}}$",
                      "$SNR_i^{\\hat{\\bf{W}},\\bf{Y}}$"))
p_snr_comp <- GLM_VB_df_w_SNR %>%
  pivot_longer(cols = c("snr_dB", "snr_dB_from_fit"), names_to = "fit", values_to = "SNR_dB") %>%
  mutate(
    fit = factor(fit, levels = c("snr_dB", "snr_dB_from_fit"), labels = c("SNR true W", "SNR est W"))
  ) %>%
  ggplot(mapping = aes(x = snr_rev_rank, y = SNR_dB, color = fit)) +
  geom_point(position = position_dodge(width = .7), size = 1.5) +
  labs(
    x = "node",
    y = TeX(r'($SNR_i$ (dB))'),
    title = TeX(r'((a) SNRs by node)'),
    subtitle = TeX(r'(from true $\bf{W}$ and estimated \hat{\bf{W}})'),
    color = ""
  ) +
  scale_color_manual(values = c("black", "firebrick1"), labels = unname(legend_labels)) +
  theme_bw() +
  theme(
    legend.position = c(.62, .85),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    plot.title=element_text(hjust=0.5),
    plot.subtitle=element_text(hjust=0.5),
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor.x = element_blank()
  ) +
  scale_x_continuous(breaks = c(1, seq(4,20,4)), minor_breaks = 1:20, limits = c(.3,20.7), expand = c(0,0))
p_snr_comp

p_ss_by_sigma <- GLM_VB_df_w_SNR %>%
  pivot_longer(cols = c("VB_fit", "GLM_fit"), names_to = "fit", values_to = "sqrd_error") %>%
  mutate(
    fit = factor(fit, levels = c("VB_fit", "GLM_fit"), labels = c("proposed VB", "naive GLM"))
  ) %>%
  ggplot(mapping = aes(x = target_sigma, y = sqrd_error, color = fit)) +
  geom_point() +
  labs(
    x = TeX(r'($\sigma_i$)'),
    y = TeX(r'($\sum_{j\ne i}(\hat{W}_{ji} - W_{ji})^2$)'),
    title = TeX(r'((b) Sum of sqrd errors)'),
    subtitle = TeX(r'(in-edge estimates by node, $\sigma_i$)'),
    color = ""
  ) +
  scale_color_manual(values = c("blue", "deeppink1")) +
  theme_bw() +
  theme(
    legend.position = c(.35, .83),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    plot.title=element_text(hjust=0.5),
    plot.subtitle=element_text(hjust=0.5),
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

p_ss_by_snr <- GLM_VB_df_w_SNR %>%
  pivot_longer(cols = c("VB_fit", "GLM_fit"), names_to = "fit", values_to = "sqrd_error") %>%
  mutate(
    fit = factor(fit, levels = c("VB_fit", "GLM_fit"), labels = c("proposed VB", "naive GLM"))
  ) %>%
  ggplot(mapping = aes(x = snr_dB, y = sqrd_error, color = fit)) +
  geom_point() +
  labs(
    x = TeX(r'($SNR_i$ (dB))'),
    y = TeX(r'($\sum_{j\ne i}(\hat{W}_{ji} - W_{ji})^2$)'),
    title = TeX(r'((b) Sum of sqrd errors)'),
    subtitle = TeX(r'(in-edge estimates by node, $SNR_i$)'),
    color = ""
  ) +
  scale_color_manual(values = c("blue", "deeppink1")) +
  theme_bw() +
  theme(
    legend.position = c(.65, .79),
    legend.title = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    plot.title=element_text(hjust=0.5),
    plot.subtitle=element_text(hjust=0.5)
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor.x = element_blank()
  )

#p_snr_comp
#p_ss_by_sigma
#p_ss_by_snr

gg <- grid.arrange(
  grobs = list(p_snr_comp + theme(plot.margin = margin(t = .1, r = .2, b = 0, l = 0, "cm")),
               p_ss_by_sigma + theme(plot.margin = margin(t = .1, r = .2, b = 0, l = .2, "cm")),
               p_ss_by_snr + theme(plot.margin = margin(t = .1, r = 1, b = 0, l = .2, "cm"))),
  nrow = 1,
  widths = c(1.2,1.1,1.2)
)

ggsave(paste0(plot_dir, "SNR_GLM.pdf"), gg, width = 8.75, height = 3.5)

p_ss_by_nsr <- GLM_VB_df_w_SNR %>%
  pivot_longer(cols = c("GLM_fit", "VB_fit"), names_to = "fit", values_to = "sqrd_error") %>%
  ggplot(mapping = aes(x = target_sigma, y = sqrd_error, color = fit)) +
  geom_point()


alpha_levels <- c(.01, .05)
TE_table_out_rows <-
  lapply(
    alpha_levels,
    function(alpha_level){
      TE_out_ALL %>%
        mutate(
          infer_edge = `p-value` < alpha_level
        ) %>%
        rowwise() %>%
        mutate(
          is_edge = !(W_true[source_node, target_node] == 0),
          correct = infer_edge == is_edge
        ) %>%
        ungroup() %>%
        group_by(T_train, is_edge) %>%
        summarize(error_rate = mean(!correct)) %>%
        ungroup() %>%
        mutate(
          rate = ifelse(is_edge, "FN", "FP"),
          names = paste0(rate, "_", T_train/1000, "k")
        ) %>%
        select(-c(is_edge,T_train,rate)) %>%
        pivot_wider(names_from = names, values_from = error_rate) %>%
        bind_cols(
          TE_out_ALL %>%
            mutate(
              infer_edge = `p-value` < alpha_level
            ) %>%
            rowwise() %>%
            mutate(
              is_edge = !(W_true[source_node, target_node] == 0),
              correct = infer_edge == is_edge
            ) %>%
            group_by(T_train) %>%
            summarize(error_rate = mean(!correct)) %>%
            ungroup() %>%
            mutate(
              names = paste0("Overall_", T_train/1000, "k")
            ) %>%
            select(-T_train) %>%
            pivot_wider(names_from = names, values_from = error_rate)
        )
    }
  ) %>%
  bind_rows %>%
  add_column(method = paste0("TE_", alpha_levels), .before = 1)


VB_table_out <- VB_out %>%
  ungroup() %>%
  mutate(
    is_edge = Wij_true != 0,
    T_train = as.numeric(run_name),
    correct = check
  ) %>%
  group_by(T_train, is_edge) %>%
  summarize(error_rate = mean(!correct)) %>%
  ungroup() %>%
  mutate(
    rate = ifelse(is_edge, "FN", "FP"),
    names = paste0(rate, "_", T_train/1000, "k")
  ) %>%
  select(-c(is_edge,T_train,rate)) %>%
  pivot_wider(names_from = names, values_from = error_rate) %>%
  bind_cols(
    VB_out %>%
      ungroup() %>%
      mutate(
        is_edge = Wij_true != 0,
        T_train = as.numeric(run_name),
        correct = check
      ) %>%
      group_by(T_train) %>%
      summarize(error_rate = mean(!correct)) %>%
      ungroup() %>%
      mutate(
        names = paste0("Overall_", T_train/1000, "k")
      ) %>%
      select(-T_train) %>%
      pivot_wider(names_from = names, values_from = error_rate)
  ) %>%
  bind_cols(
    VB_out %>%
      ungroup() %>%
      mutate(
        T_train = as.numeric(run_name),
      ) %>%
      group_by(T_train) %>%
      summarize(error_rate = sum(sqrd_error)) %>%
      ungroup() %>%
      mutate(
        names = paste0("SSE_", T_train/1000, "k")
      ) %>%
      select(-T_train) %>%
      pivot_wider(names_from = names, values_from = error_rate)
  ) %>%
  add_column(method = "proposed VB", .before = 1)


GLM_table_out <- GLM_results %>%
  ungroup() %>%
  filter(model == "stepAIC_out") %>%
  mutate(
    T_train = T_fit
  ) %>%
  group_by(T_train, is_edge) %>%
  summarize(error_rate = mean(!correct)) %>%
  ungroup() %>%
  mutate(
    rate = ifelse(is_edge, "FN", "FP"),
    names = paste0(rate, "_", T_train/1000, "k")
  ) %>%
  select(-c(is_edge,T_train,rate)) %>%
  pivot_wider(names_from = names, values_from = error_rate) %>%
  bind_cols(
    GLM_results %>%
      ungroup() %>%
      filter(model == "stepAIC_out") %>%
      mutate(
        T_train = T_fit
      ) %>%
      group_by(T_train) %>%
      summarize(error_rate = mean(!correct)) %>%
      ungroup() %>%
      mutate(
        names = paste0("Overall_", T_train/1000, "k")
      ) %>%
      select(-T_train) %>%
      pivot_wider(names_from = names, values_from = error_rate)
  ) %>%
  bind_cols(
    GLM_results %>%
      ungroup() %>%
      filter(model == "stepAIC_out") %>%
      mutate(
        T_train = T_fit,
      ) %>%
      group_by(T_train) %>%
      summarize(error_rate = sum(sqrd_error)) %>%
      ungroup() %>%
      mutate(
        names = paste0("SSE_", T_train/1000, "k")
      ) %>%
      select(-T_train) %>%
      pivot_wider(names_from = names, values_from = error_rate)
  ) %>%
  add_column(method = "naive GLM", .before = 1)


table_rows <- bind_rows(
  VB_table_out,
  TE_table_out_rows,
  GLM_table_out
)

names(table_rows)

fp_fn_rows <- table_rows[, c(1:11)] %>% as.data.frame()
rownames(fp_fn_rows) <- table_rows$method
xtable::xtable(fp_fn_rows, digits = 3)

fp_fn_rows <- table_rows[, c(12,17,13,18,14,19,15,20,16,21)] %>% as.data.frame()
rownames(fp_fn_rows) <- table_rows$method
xtable::xtable(fp_fn_rows, digits = 3)




