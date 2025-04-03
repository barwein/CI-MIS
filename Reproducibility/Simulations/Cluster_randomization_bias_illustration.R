###
# Script that illustrate the bias from using wrong adj. matrix in cluster randomization
# Cross-clusters contamination with inverse spatial distance
# Scenario III
###

# Setup -------------------------------------------------------------------


source("Reproducibility/Estimation_functions.R")

library(MisspecifiedInterference)

# Parameters and init --------------------------------------------------------------

# Clusters spec
N_clusters_fixed = 500
n_unit_each = 20
N_each_cluster_vec_fixed = rep(n_unit_each,N_clusters_fixed)

set.seed(545)
N_clusters_varied = 2000
N_each_cluster_vec_varied = sample(seq(3,7), N_clusters_varied, T)

# Sample thresholds
set.seed(586)
threshold_fixed <- runif(sum(N_each_cluster_vec_fixed))
threshold_vary <- runif(sum(N_each_cluster_vec_varied))

set.seed(5959)
base.po.noise_fixed <- runif(n = length(threshold_fixed), min = 0.5, max = 1.5)
base.po.noise_vary <- runif(n = length(threshold_vary), min = 0.5, max = 1.5)

# Sample coordinates on unit cube [0,1]^2
set.seed(5857)
x_coord_fixed = runif(N_clusters_fixed)
y_coord_fixed = runif(N_clusters_fixed)
x_coord_varied = runif(N_clusters_varied)
y_coord_varied = runif(N_clusters_varied)

clusters_coord_fixed = data.table(x_coord = x_coord_fixed, y_coord = y_coord_fixed)
clusters_coord_varied = data.table(x_coord = x_coord_varied, y_coord = y_coord_varied)

# Get distance matrix
dist_clusters_fixed = as.matrix(dist(clusters_coord_fixed, method = "manhattan"))
dist_clusters_varied = as.matrix(dist(clusters_coord_varied, method = "manhattan"))

# Theta seq
theta_seq = seq(0,0.2,0.005)
# theta_seq = seq(0,0.01,0.0005)

set.seed(2963)
Clusters_contamination_bias_fixed <- Bias_of_clusters_ATE(Nv = sum(N_each_cluster_vec_fixed),
                                                    N_clusters = N_clusters_fixed,
                                                    N_each_cluster_vec = N_each_cluster_vec_fixed,
                                                    theta_seq = theta_seq,
                                                    pz_cluster = 0.5,
                                                    threshold = threshold_fixed,
                                                    M = 10^3,
                                                    distance_matrix = dist_clusters_fixed,
                                                    base.po.noise = base.po.noise_fixed)

set.seed(2963)
Clusters_contamination_bias_varied <- Bias_of_clusters_ATE(Nv = sum(N_each_cluster_vec_varied),
                                                    N_clusters = N_clusters_varied,
                                                    N_each_cluster_vec = N_each_cluster_vec_varied,
                                                    theta_seq = theta_seq,
                                                    pz_cluster = 0.5,
                                                    threshold = threshold_vary,
                                                    M = 10^3,
                                                    distance_matrix = dist_clusters_varied,
                                                    base.po.noise = base.po.noise_vary)

























