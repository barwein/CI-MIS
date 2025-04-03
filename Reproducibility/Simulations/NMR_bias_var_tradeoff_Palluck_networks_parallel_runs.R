###
# Script that illustrate the bias-variance tradeoff of the network-misspecification-robust (NMR) estimator
# With Palluck et. al four networks as baseline
# The script is based of parallel runs which can be reporoduce in a loop
# Note that running times might be very large in a personal PC without parallelization
###

# Load libraries ----------------------------------------------------------


source("Reproducibility/Estimation_functions.R")

library(MisspecifiedInterference)

# Setup -------------------------------------------------------------------

palluck_net_list <- readRDS("Reproducibility/Data_analyses/Palluck_et_al/adj_matrix_data.RDS")
n_units <- dim(palluck_net_list[[1]])[1]

# Sample thresholds
set.seed(586)
threshold <- runif(n_units)

# Sample baseline PO
set.seed(5959)
base.po.noise <- runif(n = n_units, min = 0.5, max = 1.5)

true.adj.net <- palluck_net_list$ST.network
true.adj.name <- "ST.network"


exposures_contrast = list(c("c11","c00"),
                          c("c11","c01"),
                          c("c10","c00"),
                          c("c01","c00"),
                          c("c11","c10"))


exposures_vec = c("c11","c01","c10","c00")


# Parallel run parameter init

K <- length(palluck_net_list)

networks.combinations.mat <- matrix(data = NA,nrow=1,ncol = K)

for (k in seq(K)) {
  name.combin <- combinations(n = length(palluck_net_list),
                              r = k,
                              v = names(palluck_net_list))
  NA.mat <- matrix(NA,nrow = nrow(name.combin), ncol = K-k)
  networks.combinations.mat <- rbind(networks.combinations.mat,cbind(name.combin,NA.mat))
}
networks.combinations.mat <- networks.combinations.mat[-1,]


#######
# Replace `job_id` with values 1,...,15
# Where 15 = sum(choose(4,seq(4)))
#######


args <- commandArgs(trailingOnly = TRUE)
job_id <- NULL
if(length(args)>0){

  job_id <- as.numeric(args[1])

}else{
  stop("Error: No job id is given.")
}


cur.net.name <- networks.combinations.mat[job_id,][!is.na(networks.combinations.mat[job_id,])]
cur.A.list <- palluck_net_list[cur.net.name]

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
                    "MR_bias_var_Palluck_job_id_",
                    job_id,
                    ".csv")
write.csv(x = sim.results,
          file = file.name,
          row.names = FALSE)






