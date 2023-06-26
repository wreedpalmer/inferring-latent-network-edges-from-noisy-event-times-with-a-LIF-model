# inferring-latent-network-edges-from-noisy-event-times-with-a-LIF-model
code for dissertation chapter three and paper in progress

see paper:: [link when available]

run_simulation.Rmd creates simulated data from data model (eq. xx), contains code to create figures for publication

Running run_sim_fit.R (from run_simulation.Rmd) computes variational approximation.

Make posterior predictions for held-out individual nodes on test set with RUN_post_pred_SIM.R 

Run comparison methods of Transfer Entropy and GLM model with TE_GLM_compare.R

The other files provide supporting functions, called (sourced) in run_simulation.Rmd, RUN_post_pred_SIM.R
