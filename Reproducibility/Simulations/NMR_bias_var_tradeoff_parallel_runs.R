###
# Script that illustrate the bias-variance tradeoff of the network-misspecification-robust (NMR) estimator
# The script is based of parallel runs which can be reporoduce in a loop
# Note that running times might be very large in a personal PC without parallelization
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


# Get jaccard idx

library(kableExtra)

jacc.matrix <- matrix(NA,nrow=length(A.list),ncol=length(A.list))

for (i in seq(length(A.list))) {
  for (j in seq(length(A.list))) {
    jacc.matrix[i,j] <- jaccard_edgeset_similarity(A1 = A.list[[i]],
                                                   A2 = A.list[[j]])
  }
}

jacc.matrix[upper.tri(jacc.matrix)] <- NA
jacc.matrix <- round(jacc.matrix,3)
jacc.df <- as.data.frame(jacc.matrix)
colnames(jacc.df) <- c("$\bA^\ast$",sprintf("$\bA^%s$",seq(5)))
rownames(jacc.df) <- c("$\bA^\ast$",sprintf("$\bA^%s$",seq(5)))
options(knitr.kable.NA="")
kable(jacc.df, format = "latex", booktabs = T, digits = 3)


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


#######
# Replace `job_id` with values 1,...,63
# Where 63 = sum(choose(6,seq(6)))
#######


args <- commandArgs(trailingOnly = TRUE)
job_id <- NULL
if(length(args)>0){

  job_id <- as.numeric(args[1])

}else{
  stop("Error: No job id is given.")
}


cur.net.name <- networks.combinations.mat[job_id,][!is.na(networks.combinations.mat[job_id,])]
cur.A.list <- A.list[cur.net.name]

with.adj <- true.adj.name %in% cur.net.name

print(paste("Running Job id",job_id,"; networks used are", paste0(cur.net.name,collapse=", ")))



# Run simulations
set.seed(51552 + job_id)

sim.results <- bias_variance_full_iter_in_setup(R = 10^4,
                                                k = length(cur.net.name),
                                                # iter = 500,
                                                iter = 10^3,
                                                Pz_function = Z_ber,
                                                pz = 0.5,
                                                true.adj.mat = true.adj.net,
                                                A.list = cur.A.list,
                                                adj.mat.names = cur.net.name,
                                                with.adj = with.adj,
                                                exposures_contrast = exposures_contrast,
                                                exposures_vec = exposures_vec,
                                                threshold = threshold,
                                                base.po.noise = base.po.noise)

results.path <- "Reproducibility/Simulations/results/"
file.name <- paste0(results.path,
                    "MR_bias_var_PA_job_id_",
                    job_id,
                    ".csv")
write.csv(x = sim.results,
          file = file.name,
          row.names = FALSE)






