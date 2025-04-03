###
# Script that illustrate the bias from using wrong adj. matrix
# Here misspecification from incorrect reporting of edges (network noise)
# Scenario I
###


# Setup -------------------------------------------------------------------


source("Reproducibility/Estimation_functions.R")

library(MisspecifiedInterference)


# Parameters and init --------------------------------------------------------------

NV = 3000
PZ = 0.5
P_V = 3/NV
NET_MODELS = c("PA","ER","SW")

# Sample thresholds
set.seed(586)
threshold <- runif(NV)

# Sample baseline PO
# TODO: Check optimal DGP HERE
set.seed(5959)
base.po.noise <- runif(n = NV, min = 0.5, max = 1.5)

# Sample true networks and adj mat
set.seed(588)
true_networks <- sample_multiple_networks(net_models = NET_MODELS,
                                          Nv = NV,
                                          p_v = P_V)
# Exposure contrast
exposures_contrast = list(c("c11","c00"),
                          c("c11","c01"),
                          c("c10","c00"),
                          c("c01","c00"))


# Select true adj. matrix to run simulation on
# PA=1, ER = 2, SW = 3
net_model_name = NET_MODELS[1]
adj_true_network = true_networks$adj_mat_list[[1]]

# beta_vec is the parameters (eta) vector
beta_vec <- seq(0,0.25,0.025)

# Run simulations
set.seed(526478)
network_noise_bias_results <- Bias_of_Noisy_Network_Simulation(Nv = NV,
                                                         pz = PZ,
                                                         p_v = P_V,
                                                         net_model_name = net_model_name,
                                                         # K_vec = K_vec,
                                                         beta_vec = beta_vec,
                                                         threshold = threshold,
                                                         adj_true_network = adj_true_network,
                                                         M = 10^3,
                                                         R = 10^4,
                                                         exposures_vec = c("c11","c01","c10","c00"),
                                                         exposures_contrast = exposures_contrast,
                                                         base.po.noise = base.po.noise)



# Appendix spagehti sim ---------------------------------------------------

# Run the simulation with multiple iterations of censored networks

n.iter <- 50
# Results list
results.list <- list()

set.seed(520510)
for (n.it in seq(n.iter)) {

  cur.bias.results <- Bias_of_Noisy_Network_Simulation(Nv = NV,
                                                      pz = PZ,
                                                      p_v = P_V,
                                                      net_model_name = net_model_name,
                                                      # K_vec = K_vec,
                                                      beta_vec = beta_vec,
                                                      threshold = threshold,
                                                      adj_true_network = adj_true_network,
                                                      # M = 10^3,
                                                      M = 200,
                                                      R = 10^4,
                                                      exposures_vec = c("c11","c01","c10","c00"),
                                                      exposures_contrast = exposures_contrast,
                                                      base.po.noise = base.po.noise)

  # Add iter variable
  cur.bias.results[,net.iter := n.it]

  # update list
  results.list[[n.it]] <- cur.bias.results
}


results.dt <- rbindlist(results.list)







