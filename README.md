# Inferring latent network edges from noisy event times with a leaky integrate-and-fire (LIF) model
This is a repository for sharing R code associated with my PhD reasearch. For the full story please see chapter 3 of my dissertation "Statistical Methods for Structured Data: Analyses of Discrete Time Series and Networks," available here https://doi.org/10.7916/tfdh-kt07.

# Abstract of research
We consider the task of inferring the connections between noisy observations of events.
In our model-based approach, we consider a generative process incorporating latent dynamics that are directed by past events and the unobserved network structure. This process is based on a leaky integrate-and-fire (LIF) model from neuroscience for aggregating input and triggering events (spikes) in neural populations. Given observation data we estimate the model parameters with a novel variational Bayesian approach, specifying a highly structured and parsimonious approximation for the conditional posterior distribution of the process's latent dynamics.
This approach allows for fully interpretable inference of both the model parameters of interest and the variational parameters. Moreover, it is computationally efficient in scenarios when the observed event times are not too sparse.
We apply our methods in a simulation study and to recorded neural activity in the dorsomedial frontal cortex (DMFC) of a rhesus macaque. We assess our results based on ground truth, model diagnostics, and spike prediction for held-out nodes.

# Generative model for discrete event times

Observations from a multivariate point process of size $n$ are recorded as discretized event times in a multivariate, binary time series (matrix)
$\mathbf{Y}=(\mathrm{y}_{t,i}) \in \\{0,1\\}^{T\times n}$.
The columns of $\mathbf{Y}$ correspond to nodes in an implicit network $G = (\mathcal{V},\mathcal{E})$ with signed, weighted and directed edges $\mathcal{E}$ given by the adjacency matrix $\mathbf{W} \in \mathbb{R}^{n\times n}$.

**Connections in the latent network $G$ map past history to future observations. Specifying an explicit model for this process allows us to perform likelihood-based inference on the unobserved connections $\mathcal{E}$ from the observations $\mathbf{Y}$.**

The data generating process is a generalized state-space model for multivariate time series of binary count data. This model is adapted from the continuous-time 'leaky integrate-and-fire’ (LIF) model for neuronal dynamics, one of the simplest models describing the behavior of spiking neurons.

The adapted model incoporates variables $\mathbf{V}=(\mathrm{v}_{t,\,i}) \in \mathbb{R}^{T \times n}$, latent states corresponding to the cell membrane voltages of neurons. Observed 'event spikes' $\mathbf{y}\_t = (\mathrm{y}\_{t,1},...\,,\mathrm{y}\_{t,n})$ are Bernouilli draws, conditionally independent given  $\mathbf{v}\_t = (\mathrm{v}\_{t,1},...\,,\mathrm{v}\_{t,n})$. The latent states $\mathbf{V}$ evolve based on past observations and latent 'fluctuations' $\mathbf{Z} = (\mathrm{z}\_{t,i}) \in \mathbb{R}^{T \times n}$. The full model is given by:

![Alt text](/plots/dataModel.png?raw=true "Model")

$\mathbf{w}\_{\rightarrow i}$ is the $i\text{th}$ column of $\mathbf{W}$, giving the in-edge weights for node $i$, so that ${\mathbf{w}\_{\rightarrow i}}^\intercal \mathbf{y}\_{t} = \mathbf{W}\_{1i}\\,\mathrm{y}\_{t,1} + ... + \mathbf{W}\_{ni}\\,\mathrm{y}\_{t,n}$.

Here is a diagram of the generative model for the discrete-time observations of point process data:

![Alt text](/plots/diagram.jpg?raw=true "Model Diagram")

This figure shows the latent and observed data for a network of $n=3$ nodes over a period of length $T = 100$. It shows the model parameters for each node along with its signal-to-noise ratio (see dissertation section 3.3.9) calculated based on $\mathbf{W}$ and 100k simulated observations.

# High-level variational inference approach

With this specified model, our goal is to find values for the latent network $\mathbf{W}$ and unknown model parameters $\theta$ that maximize the marginal likelihood of the observations $\mathbf{Y}$,

$$
\mathrm{L}(\mathcal{\mathbf{W}};\mathbf{Y}) = p_{\theta, \mathbf{W}}(\mathbf{Y}) = \int p_{\theta, \mathbf{W}}(\mathbf{Y}, \mathbf{Z})\\, \mathrm{d}\mathbf{Z}.
$$

This is an intractable integral (over all the latent fluctuations $\mathbf{Z}$), so that we cannot evaluate or differentiate the marginal likelihood. 
The conditional posterior $p_{\theta, \mathbf{W}}(\mathbf{Z} | \mathbf{Y})$ is also intractable and we cannot use the EM algorithm.

Problems with similar difficulties arise often and well-known variational Bayes approaches exist. Our related approach is to use a highly structured variational approximation $q_{\phi | \mathbf{Y}} (\mathbf{Z})$ for $p_{\theta, \mathbf{W}}(\mathbf{Z} | \mathbf{Y})$, whose parameters $\phi$ we learn jointly along with the latent network $\mathbf{W}$ and unknown model parameters $\theta$ by maximizing the evidence lower bound (ELBO) of the marginal log likelihood:

![Alt text](/plots/ELBO.png?raw=true "ELBO")

By maximizing $\mathcal{L}(\theta, \mathbf{W}, \phi; \mathbf{Y})$, we minimize the Kullback–Leibler divergence between our variational approximation $q_{\phi | \mathbf{Y}} (\mathbf{Z})$ and $p_{\theta, \mathbf{W}}(\mathbf{Z} | \mathbf{Y})$.

With our proposed variational approximation $q_{\phi | \mathbf{Y}} (\mathbf{Z})$, we can take analytical gradients of a close approximation of the lower bound $\mathcal{L}(\theta, \mathbf{W}, \phi; \mathbf{Y})$ with respect to all variables $(\theta, \mathbf{W}, \phi)$. In this way we avoid noisy gradient approximations and the need for the frequently used 'reparameterization trick.'

This approach allows for fully interpretable inference of both the model parameters of interest and the variational parameters. Moreover, it is computationally efficient when the observed spike trains are not too sparse. Its complexity scales linearly with length of the observation period $T$ and number of nodes $n$, and quadratically in the length of the longest observed spike period (see section 3.3.6 for further details on complexity).

# Variational approximation

We approximate the conditional posterior distribution $p_{\theta, \mathbf{W}}(\mathbf{Z} | \mathbf{Y})$ with a parsimonious multivariate Gaussian distribution $q_{\phi | \mathbf{Y}}(\mathbf{Z})$ that imposes conditional independence and a common structure across all the inter-event periods.

Let $\mathbf{z}\_{(i,k)}$ denote the latent variables $\big(\mathrm{z}\_{t^i_k +1,i},...,\mathrm{z}\_{t^i_{k+1},i} \big)$ corresponding to the $(i,k)\text{th}$ inter-event period. 
Our approximation $q_{\phi | \mathbf{Y}}$ incorporates model parameters $(\mathbf{W}, \theta)$ and introduces variational parameters $\phi = \\{\phi_i : i = 1,...,n \\}$. It has the form

$$
q_{\phi | \mathbf{Y}}(\mathbf{Z}) = \prod\_{(i,k) \in \mathcal{X}} q_{\phi_i | \mathbf{Y}} \left(\mathbf{z}\_{(i,k)}\right ) \sim \prod\_{(i,k) \in \mathcal{X}} \mathrm{N}\Big (\mu_{\phi_i}(i,k, \mathbf{Y}),\\, \Sigma_{\phi_i}(\Delta t_k^i) \Big )
$$

where $\mathcal{X}$ is the indexed set of observed inter-event periods in $\mathbf{Y}$ with inter-event period lengths $\Delta t_k^i$. Under $q_{\phi | \mathbf{Y}}$, $\mathbf{z}\_{(i,k)}$ is normally distributed with mean $\mu_{\phi_i}(i,k, \mathbf{Y})$ and covariance $\Sigma_{\phi_i}(\Delta t_k^i)$, and $\mathbf{z}\_{(i,k)}$ is independent of $\mathbf{z}_{(i',k')}$ for all other $(i',k') \in \mathcal{X}$.

The full details of this variational approximation are in sections 3.3.2 and 3.3.3 of my dissertation. The details for computing the approximation are given in sections 3.3.4 and 3.3.5.

# Applications

I apply these methods in a simulation study and to real neural recording data. Perform the following steps in each application to compute the variational approximation:

![Alt text](/plots/computation2.png?raw=true "computation")

The details of these applications are in sections 3.4 and 3.5 of my dissertation. Here are a selected figures showing key results:

## Simulation study - comparing estimated adjacency matrix with known ground truth

![Alt text](/plots/Figure4.png?raw=true "Figure4")

## Simulation study - evaluating updates by epoch

![Alt text](/plots/Figure5.png?raw=true "Figure5")

## Simulation study - evolving estimates

![Alt text](/plots/Figure6.png?raw=true "Figure6")

## Simulation study - comparing the computed variational approximation to the posterior predictive distribution and simulated ground truth fluctuations

![Alt text](/plots/Figure11.png?raw=true "Figure11")

## Simulation study - predictive event probability estimates

![Alt text](/plots/Figure13.png?raw=true "Figure13")

## Real data application - neural recording data from dorsomedial frontal cortex (DMFC) of a rhesus macaque

![Alt text](/plots/Figure14.png?raw=true "Figure14")

## Real data application - main network inference results

![Alt text](/plots/Figure15.png?raw=true "Figure15")

# Code overview

* run_simulation.Rmd creates simulated data from the considered data model, contains code to create figures for publication

* Running run_sim_fit.R (from run_simulation.Rmd) computes the variational approximation

* Make posterior predictions for held-out individual nodes on test set with RUN_post_pred_SIM.R 

* Run comparison methods of Transfer Entropy and GLM model with TE_GLM_compare.R

* The other files provide supporting functions (are called in run_simulation.Rmd and RUN_post_pred_SIM.R)
