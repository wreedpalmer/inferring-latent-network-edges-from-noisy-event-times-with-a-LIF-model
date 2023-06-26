###### required inputs from sourcing file / environment
### ones that shouldn't change
# spike_intervals_df_restrict_nodes
# binary_spike_mat_lite_restrict_nodes
# n
### ones that may change
# T_train_cutoff
# num_epochs_run

setwd(
  ### CODE DIR HERE
)
print(paste("T_train_cutoff =", T_train_cutoff))

spike_intervals_df_train <- 
  sim_on_network$spike_intervals %>%
  filter(t_end <= T_train_cutoff)

train_indices <- spike_intervals_df_train %>% pull(index)

T_mean_inv_ob <- spike_intervals_df_train %>%
  group_by(i) %>%
  summarize(T_mean_inv = 1 / mean(T_))

# intialize variational params / sigma
phi_reparam_init <- c(
  "log_vartheta" = log(.95),
  "log_nu1" = log(.1),
  "log_nu2" = log(.25),
  "logit_u1" = gtools::logit(.9),
  "log_u2" = log(.08),  # equal to log tau
  "logit_beta1" = 0,
  "log_beta2" = log(4)
)
phi_mat_init <-
  sapply(1:n, function(node_i) get_phi_from_reparam(phi_reparam_init, ratio_bound_for_reparam = .5))

# initialize eta_vec and W matrix as zeros
#eta_init <- rep(.01, n_nodes_in_network)
#eta_init[T_mean_inv_ob$i] <- T_mean_inv_ob$T_mean_inv
eta_init <- T_mean_inv_ob$T_mean_inv
W_init <- matrix(0, nrow = (n - 1), ncol = n)

# set variables in input list to initialized values
input_list <- list(
  phi=phi_mat_init,
  eta=eta_init,
  W=W_init
)

# adam m and v variables
# start with all equal to zero
m_v_lists <- list()
m_v_lists$m <- list()
m_v_lists$v <- list()

# initialize m and v's to zero in list of lists used to keep track of them
for (cur_name in names(input_list)){
  m_v_lists$m[[cur_name]] <- input_list[[cur_name]] * 0
  m_v_lists$v[[cur_name]] <- input_list[[cur_name]] * 0
}

##### set delta and kappa (K) as (known) hyper-parameters
delta <- .975
K <- 50

#### set (tuning) parameter for gradient update steps
# moderates the step size of VI adam updates (could be set individually)
alpha_step_size = .2

# store seed for reproducibility
fitSeed <- 0
set.seed(fitSeed)

# list for storing simulation and fit
runList <- list()
runList$model_params <- list(
  n = n,
  W_true = W_true,
  delta = delta,
  eta_vec = eta_vec,
  K = K,
  sigma_vec = sigma_vec,
  num_obs = num_obs,
  sim = sim_on_network
)
runList$fit_params = list(
  spike_intervals_df_train = spike_intervals_df_train,
  binary_spike_mat_train = sim_on_network$binary_spike_train[1:T_train_cutoff, ],
  input_list_init = input_list,
  alpha_step_size = alpha_step_size,
  train_indices = train_indices,
  validation_indices = NULL,
  T_train_cutoff = T_train_cutoff,
  fitSeed = fitSeed
)

### calculate initial obj on all the data (without calculating any gradients)
#print("calculate initial obj")
#startTime <- Sys.time()
#cl <- makeCluster(13, type="FORK")
#runList$obj_init_train <- par_obj_grad_spikes(
#  cl = cl,
#  interval_indices = runList$fit_params$train_indices,
#  input_list = runList$fit_params$input_list_init,
#  spike_intervals_df = runList$fit_params$spike_intervals_df_train,
#  binary_spike_train_mat = runList$fit_params$binary_spike_mat_train,
#  delta = runList$model_params$delta,
#  K = runList$model_params$K, piecewise_approx_list=piecewise_approx_list,
#  grads = FALSE,
#  sort_on_T = FALSE
#)
#stopCluster(cl)
#print(Sys.time() - startTime)
#print(runList$obj_init_train)

### RUN adam updates
print("run adam updates")

startTime <- Sys.time()
adam_update <- update_adam_all(
  input_list = runList$fit_params$input_list_init,
  training_indices = runList$fit_params$train_indices,
  spike_intervals_df = runList$fit_params$spike_intervals_df_train,
  binary_spike_train_mat = runList$fit_params$binary_spike_mat_train,
  delta = runList$model_params$delta,
  K = runList$model_params$K, piecewise_approx_list=piecewise_approx_list,
  numNodesPar = 13,
  minibatch_count = minibatch_count,
  alpha_step_size = .2,
  m_v_lists_init = m_v_lists,
  num_epochs_run = 100,
  beta_1=0.9, beta_2=0.999, epsilon=1e-8,
  restrict_W = TRUE,
  W_lb = -.5, W_ub = .5,
  eta_lb=-.1, eta_ub=.2, transform_eta=TRUE,
  epochs_save = TRUE,
  T_train_cutoff = T_train_cutoff
)
print(Sys.time() - startTime)

runList$adam_update <- adam_update

saveRDS(
  runList,
  file=paste0(
    "xx - output/runList_sim_",
    T_train_cutoff,
    "_",
    stringr::str_extract_all(Sys.time(), "[0-9]") %>% unlist() %>% paste(collapse = ''),
    ".RData"
  )
)
