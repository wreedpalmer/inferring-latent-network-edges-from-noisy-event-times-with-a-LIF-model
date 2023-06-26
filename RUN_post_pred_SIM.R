root_dir <- "/Users/williamreedpalmer/Documents/Columbia/Research/inferring network VI/R code"
setwd(root_dir)
library(tidyverse)
library(parallel)
library(latex2exp)
library(emulator, include.only = 'quad.form')

source("polynomial_fns.R")
source("integrate_normal_fns.R")
source("get_approx_h_fns.R")
source("sim_model_fns.R")
source("comp_var_approx_fns.R")
source("samp_posterior_fns.R")
source("post_processing_fns.R")
source("post_pred_fns.R")

### LOAD runLists

runLists_filepath <- paste0(root_dir, "/0_output/sim/0_fits/")
#runList_10000 <- readRDS(paste(runLists_filepath, "runList_sim_10000_20230324022418.RData", sep = "/"))
#runList_25000 <- readRDS(paste(root_dir, "xx - output/runList_sim_25000_20230324012322.RData", sep = "/"))
runList_50000 <- readRDS(paste(runLists_filepath, "runList_sim_50000_20230414210255.RData", sep = "/"))
n <- runList_50000$model_params$n
runLists = list(#"train_10k" = runList_10000,
                #"train_25k" = runList_25000)
                "train_50k" = runList_50000)

### SET params
N_NA_sims <- 20000
N_PP_sims <- 10000
testing_length <- 2000

## run post pred pieces

# which steps to perform? and on which nodes?
run_NA_sims <- T
NA_sims_suffix <- "_Apr15_2k_20k"
run_pp_frwrd_bckwrd <- T
run_pp_suffix <- "Apr15_2k_10k"

# set runList w fit
runList <- runList_10000

# which nodes to run
held_out_nodes <- 1:20

for(runList in runLists){
  print(runList$fit_params$T_train_cutoff)
  
  # set output filepaths
  NA_sims_filepath <- paste0(root_dir, "/0_output/sim/1_R_out_NA_sims/",
                             runList$fit_params$T_train_cutoff, "/NA_sims_",
                             runList$fit_params$T_train_cutoff, NA_sims_suffix, ".RData")
  pp_out_dir <- paste0(root_dir, "/0_output/sim/3_pp/", runList$fit_params$T_train_cutoff, "/")
  
  T_train_cutoff <- runList$fit_params$T_train_cutoff
  T_test_cutoff <- T_train_cutoff + testing_length
  
  # save / load R_outs_list
  #startTime <- Sys.time()
  #s_vec_mats_list_SIM <-
  #  get_s_vec_mats_list(
  #    runList = runList,
  #    T_train_cutoff = T_train_cutoff,
  #    T_test_cutoff = T_test_cutoff,
  #    sim = TRUE,
  #    parallelize_over_nodes = TRUE,
  #    n_nodes_par = 10
  #  )
  #print(Sys.time() - startTime)
  
  # save / load R_outs_list
  if (run_NA_sims) {
    startTime <- Sys.time()
    NA_sims_list <- 
      run_NA_sims_nodes(
        runList,
        nodes_to_run = 1:20, #held_out_nodes,
        n_sims = N_NA_sims,
        T_train_cutoff = T_train_cutoff,
        T_test_cutoff = T_test_cutoff,
        parallelize_over_nodes = TRUE,
        sim = TRUE,
        n_nodes_par = 4
      )
    print(Sys.time() - startTime)
    saveRDS(NA_sims_list, NA_sims_filepath)
  } else {
    NA_sims_list <- readRDS(NA_sims_filepath)
  }
  
  
  
  
  if(run_pp_frwrd_bckwrd){
    delta <- runList$model_params$delta
    
    frwrd_bckwrd_nodes_out <-
      run_frwrd_bckwrd_nodes(
        nodes_to_run = 1:20,
        NA_sims_list = NA_sims_list,
        T_train_cutoff = T_train_cutoff,
        T_test_cutoff = T_test_cutoff,
        #s_vec_mats_list = s_vec_mats_list_SIM,
        delta = delta,
        kappa = 50,
        n_sims = N_PP_sims,
        look_ahead_limit_bckwrd = 250,
        parallelize_over_nodes = TRUE,
        n_nodes_par = 5,
        save_to_file = TRUE,
        pp_out_dir = pp_out_dir,
        file_suffix = run_pp_suffix
      )
    
  }
}

