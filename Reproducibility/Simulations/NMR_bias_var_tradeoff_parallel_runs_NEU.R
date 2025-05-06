###
# Script that compute the Number of Effective Units (NEU) of the network-misspecification-robust (NMR) estimator
###

# Load libraries ----------------------------------------------------------


source("Reproducibility/Estimation_functions.R")

library(MisspecifiedInterference)

# Setup -------------------------------------------------------------------


NV = 3000
PZ = 0.5
P_V = 3/NV

# Sample thresholds
set.seed(586)
threshold <- runif(NV)

# Sample baseline PO
set.seed(5959)
base.po.noise <- runif(n = NV, min = 0.5, max = 1.5)

# # Sample true networks and adj mat
set.seed(588)
true_networks <- sample_multiple_networks(net_models = c("PA","ER","SW"),
                                          Nv = NV,
                                          p_v = P_V)
# Name them
adj_PA <- true_networks$adj_mat_list[[1]]

# Create  noisy adj mat from true one
beta_ = 0.25
#alpha_ = beta_/50
alpha_ = beta_/100

# Sample noisy networks
set.seed(589)
noisy_adj_mat <- create_noisy_adj_mat(true_adj_mat = adj_PA,
                                      n_v = NV,
                                      alpha_ = alpha_,
                                      beta_ = beta_)

noisy_adj_mat2 <- create_noisy_adj_mat(true_adj_mat = adj_PA,
                                       n_v = NV,
                                       alpha_ = alpha_,
                                       beta_ = beta_)

noisy_adj_mat3 <- create_noisy_adj_mat(true_adj_mat = adj_PA,
                                       n_v = NV,
                                       alpha_ = alpha_,
                                       beta_ = beta_)

noisy_adj_mat4 <- create_noisy_adj_mat(true_adj_mat =adj_PA,
                                       n_v = NV,
                                       alpha_ = alpha_,
                                       beta_ = beta_)

noisy_adj_mat5 <- create_noisy_adj_mat(true_adj_mat = adj_PA,
                                       n_v = NV,
                                       alpha_ = alpha_,
                                       beta_ = beta_)


A.list <- list(adj_PA = adj_PA,
               noisy1 = noisy_adj_mat,
               noisy2 = noisy_adj_mat2,
               noisy3 = noisy_adj_mat3,
               noisy4 = noisy_adj_mat4,
               noisy5 = noisy_adj_mat5)


# Run simulation

true.adj.net <- A.list$adj_PA
true.adj.name <- names(A.list)[1]

# exposures_contrast = list(c("c11","c00"))

exposures_contrast = list(c("c11","c00"),
                          c("c11","c01"),
                          c("c10","c00"),
                          c("c01","c00"),
                          c("c11","c10"))


exposures_vec = c("c11","c01","c10","c00")


# Parallel run parameter init

K <- length(A.list)

networks.combinations.mat <- matrix(data = NA,nrow=1,ncol = K)

for (k in seq(K)) {
  name.combin <- combinations(n = length(A.list),
                              r = k,
                              v = names(A.list))
  NA.mat <- matrix(NA,nrow = nrow(name.combin), ncol = K-k)
  networks.combinations.mat <- rbind(networks.combinations.mat,cbind(name.combin,NA.mat))
}
networks.combinations.mat <- networks.combinations.mat[-1,]


# aux function for NEU -------------------------------------------------------------------


nmr_network_neu <- function(true.adj.mat, 
                            A.list, 
                            adj.mat.names, 
                            with.adj,
                            exposures_vec,
                            threshold,
                            base.po.noise, 
                            iter = 1,
                            Pz_fnction = Z_ber,
                            pz = 0.5) {
  
  # Number of units and number of networks
  n <- nrow(A.list[[1]])
  K <- length(A.list)
  
  # Initialize list to store results per iteration
  results <- vector("list", iter)
  
  for (i in seq_len(iter)) {
    # Sample a treatment assignment vector
    Z <- Pz_fnction(n, pz)
    
    # Compute exposures for each network in A.list
    expos_mat <- exposures_under_network_list(Z = Z, A.list = A.list, threshold = threshold)
    
    # For each exposure value, count units with that exposure under all networks
    NEUs <- sapply(exposures_vec, function(ck) {
      sum(specific_exposure_under_network_list(exposures_mat = expos_mat, ck = ck))
    })
    
    # Create a data.table of results for this iteration
    results[[i]] <- data.table(
      iter = i,
      K = K,
      with.true.adj = with.adj,
      exposure = exposures_vec,
      NEU = NEUs,
      adj.mat.used = paste0(adj.mat.names,collapse = ", ")
    )
  }
  
  # Combine all iterations into one data.table
  return(rbindlist(results))
}

# Run -------------------------------------------------------------------


#######
# 63 = sum(choose(6,seq(6))) --> number of networks combinations
#######


# job_id_seq <- seq(5)
job_id_seq <- seq(63)

for (job_id in job_id_seq) {
  
  # Get NEU for each combinations of networks
  
  cur.net.name <- networks.combinations.mat[job_id,][!is.na(networks.combinations.mat[job_id,])]
  cur.A.list <- A.list[cur.net.name]
  
  with.adj <- true.adj.name %in% cur.net.name
  
  print(paste("Running Job id",job_id,"; networks used are", paste0(cur.net.name,collapse=", ")))
  
  # Run simulations
  set.seed(51552 + job_id)
  
  cur.res <- nmr_network_neu(
      true.adj.mat = true.adj.net,
      A.list = cur.A.list,
      adj.mat.names = cur.net.name,
      with.adj = with.adj,
      exposures_vec = exposures_vec,
      threshold = threshold,
      base.po.noise = base.po.noise,
      # iter = 10^3,
      iter = 2,
      Pz_fnction = Z_ber,
      pz = PZ
    )
  
  # Save results
  write.csv(cur.res,
            file = paste0("Reproducibility/Simulations/neu/NMR_NEU_results_job_", job_id, ".csv"),
            row.names = FALSE)
}



