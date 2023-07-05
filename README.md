# Inferring latent network edges from noisy event times with a leaky integrate-and-fire (LIF) model
Code for PhD reasearch (see chapter 3 of dissertation here https://doi.org/10.7916/tfdh-kt07)

# Abstract of research
We consider the task of inferring the connections between noisy observations of events.
In our model-based approach, we consider a generative process incorporating latent dynamics that are directed by past events and the unobserved network structure. This process is based on a leaky integrate-and-fire (LIF) model from neuroscience for aggregating input and triggering events (spikes) in neural populations. Given observation data we estimate the model parameters with a novel variational Bayesian approach, specifying a highly structured and parsimonious approximation for the conditional posterior distribution of the process's latent dynamics.
This approach allows for fully interpretable inference of both the model parameters of interest and the variational parameters. Moreover, it is computationally efficient in scenarios when the observed event times are not too sparse.
We apply our methods in a simulation study and to recorded neural activity in the dorsomedial frontal cortex (DMFC) of a rhesus macaque. We assess our results based on ground truth, model diagnostics, and spike prediction for held-out nodes.

$$
x \mapsto f(x)
$$

![Alt text](/plots/diagram.jpg?raw=true "Optional Title")

run_simulation.Rmd creates simulated data from data model (eq. xx), contains code to create figures for publication

Running run_sim_fit.R (from run_simulation.Rmd) computes variational approximation.

Make posterior predictions for held-out individual nodes on test set with RUN_post_pred_SIM.R 

Run comparison methods of Transfer Entropy and GLM model with TE_GLM_compare.R

The other files provide supporting functions, called (sourced) in run_simulation.Rmd, RUN_post_pred_SIM.R
