
library(igraph)
library(foreach)
library(doParallel)
library(data.table)
library(mgcv)
library(parallel)
library(gtools)

# install.packages("devtools")
# TEST
#devtools::install_github("barwein/Misspecified_Interference")
library(MisspecifiedInterference)

# Simulations -----------------------------------------------------------

generate_exposures_threshold <- function(A, Z, threshold,deg=NULL){
  # Compute num neighbors
  if(is.null(deg)){
    deg <- apply(A, 1, sum)
  }
  # Compute signs
  indirect_sign <- as.numeric((A %*% Z) > threshold*deg)

  # Save exposures
  exposures <- vector(length = nrow(A))

  exposures[Z*indirect_sign==1] <- "c11"
  exposures[(1-Z)*indirect_sign==1] <- "c01"
  exposures[Z*(1-indirect_sign)==1] <- "c10"
  exposures[(1-Z)*(1-indirect_sign)==1] <- "c00"

  return(exposures)
}

generate_po <- function(exposures, base.noise = NULL){

  if(is.null(base.noise)){
    n <- length(exposures)
    y.00 <- runif(n = n, min = 0.5, max = 1.5)
  } else{
    y.00 <- base.noise
  }
  # Y <- 2*y.00*(exposures=="c11") +
  #   1.25*y.00*(exposures=="c01") +
  #   1.5*y.00*(exposures=="c10") +
  #   1*y.00*(exposures=="c00")

  Y <- (y.00+1)*(exposures=="c11") +
    (y.00+0.25)*(exposures=="c01") +
    (y.00+0.5)*(exposures=="c10") +
    y.00*(exposures=="c00")

  return(Y)
}


Z_ber <- function(n,p){
  rbinom(n = n, size = 1, prob = p)
}

Z_ber_clusters <- function(N_clusters,
                           N_each_cluster_vec,
                           p){

  treated_cluster_indicator = rep(0, N_clusters)
  treated_cluster_idx = sample(x = seq(N_clusters),
                               size = N_clusters %/% (1/p),
                               replace = FALSE)
  treated_cluster_indicator[treated_cluster_idx] <- 1
  unit_level_treatment_vec = rep(treated_cluster_indicator, N_each_cluster_vec)
  return(unit_level_treatment_vec)
}

compute_exposure_prob <- function(Pz_function, pz_func_args, R, A, threshold){

  # Sample R Z-vector and save in a matrix (nv*R)
  Z_mat <- replicate(R, do.call(Pz_function, pz_func_args))

  # Convert treatments to exposures
  expos_mat <- apply(Z_mat, 2, function(x){
    generate_exposures_threshold(A = A,Z = x, threshold = threshold)
  })

  # Compute prob. using additive smoothing
  p_11 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c11")+1)/(R+1)})
  p_01 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c01")+1)/(R+1)})
  p_10 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c10")+1)/(R+1)})
  p_00 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c00")+1)/(R+1)})

  return(list("p11" = p_11, "p01" = p_01,
              "p10" = p_10, "p00" = p_00))

}

binom_prob <- function(a,b,p,n){
  return(pbinom(q = b,size = n, prob = p)-pbinom(q = a,size = n, prob = p))
}

compute_exposure_prob_analytical_ber <- function(A, p, threshold){

  Nv <- nrow(A)

  # Compute degrees vectors
  deg_vec <- apply(A, 2, sum)

  # Compute general threshold and round below
  general_threshold <- floor(deg_vec*threshold)

  # Compute exposure prob.

  p_11 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){p*binom_prob(a=x[2], b=x[1], p=p, n=x[1])})

  p_01 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){(1-p)*binom_prob(a=x[2], b=x[1], p=p, n=x[1])})

  p_10 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){p*binom_prob(a=-1, b=x[2], p=p, n=x[1])})

  p_00 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){(1-p)*binom_prob(a=-1, b=x[2], p=p, n=x[1])})

  # Return results
  return(list("p11" = p_11, "p01" = p_01,
              "p10" = p_10, "p00" = p_00))

}

compute_joint_exposure_prob <- function(Pz_function, nv, R, A.true, A.sp, p, threshold, exposures_vec){
  
  deg_A.true <- apply(A.true, 1 , sum)
  deg_A.sp <- apply(A.sp, 1 , sum)
  
  n.exp <- length(exposures_vec)
  
  joint.count.mat <- as.data.frame(matrix(0,nrow = nv, ncol = n.exp^2))
  
  all.exposure.combination <- rbind(rep(exposures_vec,each=4),rep(exposures_vec,times=4))
  names(joint.count.mat) <- apply(all.exposure.combination, 2, function(x){paste0(x[1],",",x[2])})
  # Exposure format is (true.exp, sp.exp)
  
  for (m in seq(R)) {
    # sample treatment vector
    z <- Pz_function(n=nv,p=p)
    # Compute exposures
    exposures_A.true <- generate_exposures_threshold(A = A.true,
                                                     Z = z,
                                                     threshold = threshold,
                                                     deg = deg_A.true)
    exposures_A.sp <- generate_exposures_threshold(A = A.sp,
                                                   Z = z, 
                                                   threshold = threshold,
                                                   deg = deg_A.sp)
    # update joint exposures counter 
    for (i in seq(nv)){
      joint.count.mat[i,paste0(exposures_A.true[i],",",exposures_A.sp[i])] <- 
        joint.count.mat[i,paste0(exposures_A.true[i],",",exposures_A.sp[i])] + 1
      
    }
  }
  
  # Return prob.mat (Nv*(n.exposures^2) dimension)
  return((joint.count.mat + 1) / (R+1))
}


compute_true_exposure_prob <- function(network, p_z){

  deg_vec <- degree(network)

  true_p11 <- p_z*(1-((1-p_z)^deg_vec))
  true_p01 <- (1-p_z)*(1-((1-p_z)^deg_vec))
  true_p10 <- p_z*((1-p_z)^deg_vec)
  true_p00 <- (1-p_z)*((1-p_z)^deg_vec)

  return(list(true_p11 = true_p11,
              true_p01 = true_p01,
              true_p10 = true_p10,
              true_p00 = true_p00))
}

compute_exposure_misclassification <- function(true_exposure,
                                               noisy_exposure){
  m.c11 <- sum(noisy_exposure == "c11" & true_exposure != "c11")
  m.c01 <- sum(noisy_exposure == "c01" & true_exposure != "c01")
  m.c10 <- sum(noisy_exposure == "c10" & true_exposure != "c10")
  m.c00 <- sum(noisy_exposure == "c00" & true_exposure != "c00")
  return(c(m.c11 = m.c11,
           m.c01 = m.c01,
           m.c10 = m.c10,
           m.c00 = m.c00))
}

create_noisy_adj_mat <- function(true_adj_mat,
                                 n_v,
                                 alpha_,
                                 beta_){

  # Save only lower triangular matrix

  lower_tri_idx <- lower.tri(true_adj_mat)

  noisy_adj_mat <- true_adj_mat

  noisy_adj_mat[!lower_tri_idx] <- NA

  # Remove existing ties with prob. beta_
  # Add ties with prob. alpha_

  ind_ties <- which(noisy_adj_mat==1)
  ind_no_ties <- which(noisy_adj_mat==0)
  noisy_adj_mat[ind_ties] <- rbinom(length(ind_ties),1,1-beta_)
  noisy_adj_mat[ind_no_ties] <- rbinom(length(ind_no_ties),1,alpha_)

  # Compute the full noisy matrix from the lower tri one

  # diag is 0
  diag(noisy_adj_mat) <- 0

  # Force the upper tri be the same as the lower tri

  noisy_adj_mat <- as.matrix(Matrix::forceSymmetric(noisy_adj_mat,uplo="L"))

  return(noisy_adj_mat)
}


sample_network <- function(net_mod,
                           Nv,
                           p_v = NULL){
  if(net_mod == "PA"){ # Preferential attachment
    return(sample_pa(n = Nv, power = 1, directed = FALSE))
  }
  if(net_mod == "ER"){ # Erdős–Rényi random graph
    return(sample_gnp(n = Nv, p = p_v))
  }
  if(net_mod == "SW"){ # Small-world graph
    G <- sample_smallworld(dim = 1,
                           size = Nv,
                           nei = 2,
                           p = 0.1)
    return(simplify(G, remove.loops = TRUE))
  }
}

sample_multiple_networks <- function(net_models,
                                     Nv,
                                     p_v = NULL){
  network_list <- list()
  adj_mat_list <- list()
  for (n_net in seq(length(net_models))){
    cur_net <- sample_network(net_mod = net_models[n_net],
                              Nv = Nv,
                              p_v = p_v)
    cur_adj_mat <- as.matrix(as_adjacency_matrix(cur_net))

    network_list[[n_net]] <- cur_net
    adj_mat_list[[n_net]] <- cur_adj_mat
  }
  return(list(network_list = network_list,
              adj_mat_list = adj_mat_list))
}



n_joint_exposures <- function(expos_A, expos_obs){
  n.11 <- sum(expos_A=="c11" & expos_obs=="c11")
  n.01 <- sum(expos_A=="c01" & expos_obs=="c01")
  n.10 <- sum(expos_A=="c10" & expos_obs=="c10")
  n.00 <- sum(expos_A=="c00" & expos_obs=="c00")
  return(c(n.11, n.01, n.10, n.00))
}


exposures_under_network_list <- function(Z, A.list, threshold){
  # for m adj. matrices in A.list and n units, this function compute the nXm matrix of exposures under
  # each adj. matrix for treatment vector Z
  return(
    sapply(A.list, function(A){generate_exposures_threshold(A=A,Z=Z,threshold=threshold)})
  )}


specific_exposure_under_network_list <- function(exposures_mat, ck){
  # for some exposures matrix, this function return a binary (1/0) vector
  # of units that have exposure ck under all adj. matrics (rows in exposure_mat)
  return(
    apply(exposures_mat,1,function(x){as.numeric(all(x==ck))})
  )}


compute_prob_indicator <- function(R, n, Pz_function, pz_func_args, A.list, exposures_vec, threshold){
  # A function that returns a list of nXr indicator matrices (one for each exposures),
  # that indicate if unit i in sample r is exposed to some ck under all A \in A.list
  # The output is used for computation of exposures probabilities
  
  
  # Init a list of nXR indicator matrices 
  Ind.mat.list <- vector("list",length(exposures_vec))
  for (i in seq(length(exposures_vec))) {
    assign(paste0("I.",exposures_vec[i]), matrix(NA, nrow = n, ncol = R))
    Ind.mat.list[[i]] <- get(paste0("I.",exposures_vec[i]))
  }
  names(Ind.mat.list) <- exposures_vec
  
  # Sample Z ~ PZ R times and update the indicators matrices
  for (r in seq(R)) {
    # cur.Z <- Pz_function(n=n, p=pz)
    cur.Z <- do.call(Pz_function, pz_func_args)
    cur.exposure.mat <- exposures_under_network_list(Z=cur.Z,
                                                     A.list=A.list,
                                                     threshold=threshold)
    for (ck in exposures_vec) {
      # Update indic?tor matrices for each exposures and r
      Ind.mat.list[[ck]][,r] <- specific_exposure_under_network_list(exposures_mat=cur.exposure.mat,
                                                                     ck=ck)
    }
  }
  return(Ind.mat.list)
}


compute_prob_matrices <- function(ind.mat.list, exposures_contrast, exposures_vec){
  
  # From the list of indicator matrices, vector of exposures and relevant contrast,
  # the function return a list of all relevant probability matrices
  
  n <- nrow(ind.mat.list[[1]])
  R <- ncol(ind.mat.list[[1]])
  # Init list that will contain all relevant probability matrices
  P_list <- vector(mode = "list", length = length(exposures_vec) + length(exposures_contrast))
  
  # Update P.k matrices
  for (i in seq(length(exposures_vec))) {
    ind.mat <- ind.mat.list[[exposures_vec[i]]]
    # Use additive smoothing for the diagonal values
    P_list[[i]] <- (ind.mat %*% t(ind.mat) + diag(1,n)) / (R+1) 
  }
  
  if(length(exposures_contrast) > 0){
    # Update P.kl matrices
    for (j in seq(length(exposures_contrast))) {
      ck <- exposures_contrast[[j]][1]
      cl <- exposures_contrast[[j]][2]
      ind.mat.k <- ind.mat.list[[ck]]
      ind.mat.l <- ind.mat.list[[cl]]
      P_list[[length(exposures_vec)+j]] <- (ind.mat.k %*% t(ind.mat.l)) / R
    }
  }
  names(P_list) = c(exposures_vec,
                    sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])}))
  
  return(P_list)
}


Get_prob_matrices_list <- function(R, n,
                                   Pz_function,
                                   pz_func_args,
                                   A.list,
                                   exposures_contrast, exposures_vec,
                                   threshold,
                                   Palluck.et.al = FALSE,
                                   schid = NULL){
  # Wrapper function that take as input all the adj. matrix list, R (number of iteration),
  # exposures, exposures contrasts list, and Z sampler, 
  # and comptute the list of probability matrices for the later use by MR estimator
  
  # First, get the indicator matrix list (list with length equal to exposures_vec,
  #                                         one matrix for each exposure)
  if(Palluck.et.al){
    ind.mat.list <- compute_prob_indicator_palluck(R = R, n = n,
                                                   Pz_function = Pz_function,
                                                   pz_func_args = pz_func_args,
                                                   # schid = schid,
                                                   A.list = A.list,
                                                   exposures_vec = exposures_vec)  
  } else{
    ind.mat.list <- compute_prob_indicator(R = R, n = n,
                                           Pz_function = Pz_function,
                                           pz_func_args = pz_func_args,
                                           A.list = A.list,
                                           exposures_vec = exposures_vec,
                                           threshold = threshold)
  }
  
  # Compute the relevant prob. matrices
  Prob_matrices_list <- compute_prob_matrices(ind.mat.list = ind.mat.list,
                                              exposures_contrast = exposures_contrast,
                                              exposures_vec = exposures_vec)
  
  return(Prob_matrices_list)
}



compute_MR_PO_estimator <- function(Y.obs, expos.obs.mat, ck, P.k){
  # Function that estimate the mean potential outcome using the MR estimator
  
  n <- length(Y.obs)
  # Get indicator of exposure ck vector
  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)
  
  # Get probabilities from the diagonal of P.k
  Pii <- diag(P.k)
  
  # Estimate using both HT and Hajek
  
  HT_esti <- (n^{-1})*(sum(ind.k*Y.obs/Pii))
  
  hajek_esti <- (sum(ind.k*Y.obs/Pii)) / (sum(ind.k/Pii))
  
  naive_esti <- sum(ind.k*Y.obs) / sum(ind.k)
  
  return(list(ht_esti = HT_esti, hajek_esti = hajek_esti, naive_esti = naive_esti))
}

compute_var_MR_PO_estimator <- function(Y.obs, expos.obs.mat, ck, P.k, estimate.n = FALSE){
  # Function that compute the variance of the mean potential outcome  MR estimator
  
  n <- length(Y.obs)
  # Get indicator of exposure ck vector
  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)
  
  # Get probabilities from the diagonal of P.k
  Pii <- diag(P.k)
  
  # Whether to estimate n or use given n (HT vs Hajek)
  n.hat <- ifelse(estimate.n,sum(ind.k/Pii),n)
  
  # First sum of the formula
  sum.1 <- (n.hat^{-2})*(sum(ind.k*
                               (1-Pii)*
                               ((Y.obs/Pii)^2)
  ))  
  
  sum.2 = 0; sum.3 = 0
  
  # Iteratively update second and third of the formula
  for (i in seq(n)) {
    for (j in setdiff(seq(n),i)){
      if(P.k[i,j] != 0 & ind.k[i]*ind.k[j]==1){
        # if(ind.k[i]*ind.k[j]==1){
        cur.sum = ((P.k[i,j] - Pii[i]*Pii[j])/(P.k[i,j]))*
          ((Y.obs[i]*Y.obs[j])/(Pii[i]*Pii[j]))
        sum.2 = sum.2 + cur.sum
      }
      if(P.k[i,j] == 0){
        cur.sum = (ind.k[i]*(Y.obs[i]^2)/(2*Pii[i])) +
          (ind.k[j]*(Y.obs[j]^2)/(2*Pii[j]))
        sum.3 = sum.3 + cur.sum
      }
    }  
  }
  # Scale the sums
  sum.2 = (n.hat^{-2})*sum.2
  sum.3 = (n.hat^{-2})*sum.3
  # Return variance estimator
  return(sum.1 + sum.2 + sum.3)
}


compute_cov_MR_PO_estimator <- function(Y.obs, expos.obs.mat, ck, cl,
                                        P.k, P.l, P.kl,
                                        estimate.n = FALSE){
  
  n <- length(Y.obs)
  # Get indicator of exposure ck vector
  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)
  ind.l <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = cl)
  
  # Whether to estimate n or use given n (HT vs Hajek)
  # n.hat <- ifelse(estimate.n,sum(ind.k/Pii),n)
  n.hat <- n
  
  sum.1 = 0; sum.2 = 0
  
  for (i in seq(n)) {
    for (j in seq(n)) {
      if(P.kl[i,j] > 0 & ind.k[i]*ind.l[j] == 1 & i != j){
        cur.sum = ((P.kl[i,j] - P.k[i,i]*P.l[j,j])/(P.kl[i,j]))*
          ((Y.obs[i]*Y.obs[j])/(P.k[i,i]*P.l[j,j]))
        sum.1 = sum.1 + cur.sum
      }
      if(P.kl[i,j] == 0){
        cur.sum = (ind.k[i]*(Y.obs[i]^2)/(2*P.k[i,i])) +
          (ind.l[j]*(Y.obs[j]^2)/(2*P.l[j,j])) 
        sum.2 = sum.2 + cur.sum
      }
    }
  }
  sum.1 = (n.hat^{-2})*sum.1
  sum.2 = (n.hat^{-2})*sum.2
  
  return(sum.1 - sum.2)
}

# NOTE that each prob. matrix only have to be computed ONCE! that will reduce running time!


Hajek_residualized_Y <- function(Y.obs, ck, cl, expos.obs.mat, hajek_esti.ck, hajek_esti.cl){
  # Get linear residualized values from Hajek estimator 
  # only for units with exposure ck under all adj. matrix list.
  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)
  ind.l <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = cl)
  y.resid <- ind.k*(Y.obs - hajek_esti.ck) + ind.l*(Y.obs - hajek_esti.cl) 
  return(y.resid)
}


MR_CE_wrapper <- function(Y.obs, expos.obs.mat, ck, cl,
                          P.k, P.l, P.kl, 
                          estimate.n = FALSE,
                          compute.variance = TRUE){
  # Wrapper function that compute point estimate and variance of causal effect \tau(ck,cl)
  
  # Get point estimates using HT and Hajek
  mu.ck <- compute_MR_PO_estimator(Y.obs = Y.obs,
                                   expos.obs.mat = expos.obs.mat,
                                   ck = ck,
                                   P.k = P.k)
  
  mu.cl <- compute_MR_PO_estimator(Y.obs = Y.obs,
                                   expos.obs.mat = expos.obs.mat,
                                   ck = cl,
                                   P.k = P.l)
  
  ht_ce <- mu.ck$ht_esti - mu.cl$ht_esti
  hajek_ce <- mu.ck$hajek_esti - mu.cl$hajek_esti  
  naive_diff_ce <- mu.ck$naive_esti - mu.cl$naive_esti
  
  if(compute.variance){
    # Get variance estimator for HT
    var.ht.ck <- compute_var_MR_PO_estimator(Y.obs = Y.obs,
                                             expos.obs.mat = expos.obs.mat,
                                             ck = ck,
                                             P.k = P.k,
                                             estimate.n = estimate.n)
    
    var.ht.cl <- compute_var_MR_PO_estimator(Y.obs = Y.obs,
                                             expos.obs.mat = expos.obs.mat,
                                             ck = cl,
                                             P.k = P.l,
                                             estimate.n = estimate.n)
    
    cov.ht.ckcl <- compute_cov_MR_PO_estimator(Y.obs = Y.obs,
                                               expos.obs.mat = expos.obs.mat,
                                               ck = ck,
                                               cl = cl, 
                                               P.k = P.k,
                                               P.l = P.l,
                                               P.kl = P.kl,
                                               estimate.n = estimate.n)
    
    var_ht_ce <- var.ht.ck + var.ht.ck - 2*cov.ht.ckcl
    
    # Get variance estimator for Hajek
    # First need to residualize Y.obs and then use mentioned variance estimator
    Y.resid <- Hajek_residualized_Y(Y.obs = Y.obs,
                                    ck = ck,
                                    cl = cl, 
                                    expos.obs.mat = expos.obs.mat,
                                    hajek_esti.ck = mu.ck$hajek_esti,
                                    hajek_esti.cl = mu.cl$hajek_esti)
    
    var.hajek.ck <- compute_var_MR_PO_estimator(Y.obs = Y.resid,
                                                expos.obs.mat = expos.obs.mat,
                                                ck = ck,
                                                P.k = P.k,
                                                estimate.n = estimate.n)
    
    var.hajek.cl <- compute_var_MR_PO_estimator(Y.obs = Y.resid,
                                                expos.obs.mat = expos.obs.mat,
                                                ck = cl,
                                                P.k = P.l,
                                                estimate.n = estimate.n)
    
    cov.hajek.ckcl <- compute_cov_MR_PO_estimator(Y.obs = Y.resid,
                                                  expos.obs.mat = expos.obs.mat,
                                                  ck = ck,
                                                  cl = cl, 
                                                  P.k = P.k,
                                                  P.l = P.l,
                                                  P.kl = P.kl,
                                                  estimate.n = estimate.n)
    
    var_hajek_ce <- var.hajek.ck + var.hajek.cl - 2*cov.hajek.ckcl
  }
  if(!compute.variance){
    var_ht_ce = var_hajek_ce = NA
  }
  # Return results
  return(list(ht_ce = ht_ce, hajek_ce = hajek_ce,
              naive_diff_ce = naive_diff_ce, 
              var_ht_ce = var_ht_ce, var_hajek_ce = var_hajek_ce))
}

MR_CE_estimator <- function(Z.obs, Y.obs, A.list,
                            exposures_contrast, exposures_vec,
                            Prob_matrices_list,
                            threshold,
                            estimate.n = FALSE,
                            compute.variance = TRUE){
  
  expos.obs.mat <- exposures_under_network_list(Z = Z.obs, A.list = A.list, threshold = threshold)
  
  ce.results <- vector("list", length(exposures_contrast))
  names(ce.results) <- sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])})
  
  for (expos in exposures_contrast) {
    ck <- expos[1]
    cl <- expos[2]
    cntrast_name <- paste0(ck,"-",cl)
    P.k <- Prob_matrices_list[[ck]]
    P.l <- Prob_matrices_list[[cl]]
    P.kl <- Prob_matrices_list[[cntrast_name]]
    ce_estimate <- MR_CE_wrapper(Y.obs = Y.obs,
                                 expos.obs.mat = expos.obs.mat,
                                 ck = ck,
                                 cl = cl,
                                 P.k = P.k,
                                 P.l = P.l,
                                 P.kl = P.kl,
                                 estimate.n = estimate.n,
                                 compute.variance = compute.variance)
    ce.results[[cntrast_name]] <- ce_estimate
  }
  return(ce.results)
}



mod_rbind <- function(dt1, dt2){
  rbindlist(list(dt1,dt2))
}


Get_potential_outcomes_matrix <- function(threshold, base.noise=NULL){
  if(is.null(base.noise)){
    n <- length(exposures)
    y.00 <- runif(n = n, min = 0.5, max = 1.5)
  } else{
    y.00 <- base.noise
  }
  
  # Y <- as.data.frame(cbind(2*y.00, 1.25*y.00, 1.5*y.00, 1*y.00))
  Y <- as.data.frame(cbind(y.00 + 1, y.00 + 0.25, y.00 + 0.5, y.00))
  
  names(Y) <- c("c11","c01","c10","c00")
  
  return(Y)
}


Get_bias_weight <- function(sp.exposure,
                            A.sp.prob.vec,
                            Joint.prob.mat,
                            exposures_vec){
  
  relevant.cols <- paste0(exposures_vec,",",sp.exposure)
  which.same.exposure <- which(exposures_vec==sp.exposure)
  # Compute weight
  q.mat <- Joint.prob.mat[,relevant.cols] / A.sp.prob.vec
  # Adjust for j=k
  q.mat[,which.same.exposure] <- q.mat[,which.same.exposure] - 1
  
  names(q.mat) <- exposures_vec
  
  return(q.mat)
}

Get_bias_bounds_probs <- function(noisy.prob.mat.list,
                                   Joint.prob.mat,
                                   exposures_vec){
  cond.prob.list <- list("c11"=NA, "c01"=NA, "c10"=NA, "c00"=NA)
  for (expos in exposures_vec) {
    name.pasted <- paste0(expos,",",expos)
    joint.prob <- Joint.prob.mat[,name.pasted]
    marginal.sp.prob <- diag(noisy.prob.mat.list[[expos]])
    cond.prob <- joint.prob/marginal.sp.prob
    # Update list
    cond.prob.list[expos] <- sum(1-cond.prob)  
  }
  return(cond.prob.list)
}

Get_bias_exact_and_bound <- function(noisy.prob.mat.list,
                                     true_adj_mat,
                                     noisy_adj_mat,
                                     exposures_vec,
                                     exposures_contrast, 
                                     threshold,
                                     base.po.noise,
                                     Nv,
                                     R,
                                     pz,
                                     kapa){
  # Function that compute the HT estimator the exact bias and the bounds
  
  # Get joint exposures probabilities matrix  (A.true.exposure, A.sp.exposure)
  Joint.prob.mat <- compute_joint_exposure_prob(Pz_function = Z_ber,
                                                nv = Nv,
                                                R = R,
                                                A.true = true_adj_mat, 
                                                A.sp = noisy_adj_mat,
                                                p = pz,
                                                threshold = threshold,
                                                exposures_vec = exposures_vec)
  # Get PO matrix
  Y.mat <- Get_potential_outcomes_matrix(threshold = threshold,
                                         base.noise = base.po.noise)
  
  # Get bounds weights
  bounds.weight.list <- Get_bias_bounds_probs(noisy.prob.mat.list = noisy.prob.mat.list,
                                              Joint.prob.mat = Joint.prob.mat,
                                              exposures_vec = exposures_vec)
  
  n.contrast <- length(exposures_contrast)
  # Results matrix
  bias.mat <- data.frame(ce_contrast = sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])}),
                         exact.bias = rep(NA, n.contrast),
                         maximal.bias = rep(NA, n.contrast))
  # Compute exact bias for each exposures contrast
  for (nc in seq(n.contrast)){
    cntrst <- exposures_contrast[[nc]]
    expsr.l <- cntrst[1]; expsr.r <- cntrst[2];
    
    q.mat.l <- Get_bias_weight(sp.exposure = expsr.l,
                               A.sp.prob.vec = diag(noisy.prob.mat.list[[expsr.l]]),
                               Joint.prob.mat = Joint.prob.mat,
                               exposures_vec = exposures_vec)
    
    q.mat.r <- Get_bias_weight(sp.exposure = expsr.r,
                               A.sp.prob.vec = diag(noisy.prob.mat.list[[expsr.r]]),
                               Joint.prob.mat = Joint.prob.mat,
                               exposures_vec = exposures_vec)
    
    q.diff <- q.mat.l - q.mat.r
    
    # exact bias
    bias.mat[nc,2] <- (Nv^{-1})*sum(q.diff*Y.mat)
    # bounds of bias
    bias.mat[nc,3] <- (2*kapa/Nv)*(bounds.weight.list[[expsr.l]]+bounds.weight.list[[expsr.r]])
    
  }
  
  # Return Bias(c_l,c_k; A^sp) and bounds of all causal contrasts
  return(bias.mat)
}


One_network_CE_estimation <- function(true_adj_mat,
                                      noisy_adj_mat,
                                      network_model,
                                      M,
                                      Nv,
                                      pz,
                                      R,
                                      param,
                                      # alpha_,
                                      # beta_,
                                      threshold,
                                      exposures_vec,
                                      exposures_contrast,
                                      base.po.noise = NULL,
                                      jaccard_idx,
                                      compute.exact.bias = TRUE,
                                      kapa = 1.5+1){
  # Get prob. matrices for current censored network
  prob.mat <- Get_prob_matrices_list(R = R, n = Nv,
                                     Pz_function = Z_ber, pz_func_args = list(n=Nv,p=pz),
                                     A.list = list(noisy_adj_mat),
                                     exposures_contrast = exposures_contrast,
                                     exposures_vec = exposures_vec,
                                     # exposure_func = generate_exposures_threshold,
                                     threshold = threshold
                                     # exposure_func_args = list(threshold)
                                     )
  # Estimate CE M times
  sim_results <- mclapply(seq(M), function(m){
    # Sample treatment vector
    Z <- Z_ber(n = Nv, p = pz)
    # Compute true exposures
    true_expos <- generate_exposures_threshold(A = true_adj_mat,
                                               Z = Z,
                                               threshold = threshold)
    # Compute observed outcomes (using TRUE adj mat)
    true_outcomes <- generate_po(true_expos, base.po.noise)
    # Compute exposures
    noisy_expos <- generate_exposures_threshold(A = noisy_adj_mat,
                                                Z = Z,
                                                threshold = threshold)

    CE_estimate <- rbindlist(MR_CE_estimator(Z.obs = Z,
                                             Y.obs = true_outcomes,
                                             A.list = list(noisy_adj_mat),
                                             exposures_contrast = exposures_contrast,
                                             exposures_vec = exposures_vec,
                                             Prob_matrices_list = prob.mat,
                                             exposure_func = generate_exposures_threshold,
                                             exposure_func_args = list(threshold),
                                             compute.variance = FALSE))[,1:2]

    # Compute number of exposure mis-classifications
    expos_misclass <- compute_exposure_misclassification(true_exposure = true_expos,
                                                         noisy_exposure = noisy_expos)

    # Update results DT
    CE_estimate[,`:=`(ce_contrast = sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])}),
                      network_model = network_model,
                      param = paste0(param, collapse = ", "),
                      iter = m,
                      n.expos.misclass = paste0(names(expos_misclass),"=",expos_misclass, collapse = "; "),
                      jac = jaccard_idx
    )]

    CE_estimate
  },
  mc.cores = 1)
  
  sim.results.binded <- rbindlist(sim_results)
  
  if(compute.exact.bias){
    exact.bias.mat <- Get_exact_bias(noisy.prob.mat.list = prob.mat,
                                     true_adj_mat = true_adj_mat,
                                     noisy_adj_mat = noisy_adj_mat,
                                     exposures_vec = exposures_vec,
                                     exposures_contrast = exposures_contrast,
                                     threshold = threshold,
                                     base.po.noise = base.po.noise,
                                     Nv = Nv,
                                     R = R, pz = pz,
                                     kapa = kapa)
    sim.results.binded[,exact.bias := rep(c(exact.bias.mat[,2]),times=M)]
    sim.results.binded[,bounds.bias := rep(c(exact.bias.mat[,3]),times=M)]
  }
  
  return(sim.results.binded)
}

Cluster_randomization_CE_estimation <- function(N_clusters,
                                                N_each_cluster_vec,
                                                theta_,
                                                Nv,
                                                pz_cluster,
                                                # Z,
                                                threshold,
                                                M,
                                                distance_matrix=NULL,
                                                base.po.noise = NULL){

  # Compute exposures assuming clusters network structure
  # cluster_expos <- generate_exposures_threshold(A = clusters_adj_mat,
  #                                               Z = Z)

  # Get SBM Q matrix
  if(is.null(distance_matrix)){
    cur_Q <- generate_Q_matrix(N_clusters = N_clusters,
                               within_prob_vec = rep(1, N_clusters),
                               between_prob_vec = rep(theta_, choose(N_clusters,2)))
  } else{
    # Get lower tri
    dist_vec = distance_matrix[lower.tri(distance_matrix)]
    cur_Q <- generate_Q_matrix(N_clusters = N_clusters,
                               within_prob_vec = rep(1, N_clusters),
                               between_prob_vec = exp(-dist_vec/(mean(dist_vec)*theta_)))
  }

  # Sample SBM
  cur_sbm <- sample_sbm(n = Nv,
                        pref.matrix = cur_Q,
                        block.sizes = N_each_cluster_vec)

  sbm_adj_mat <- as.matrix(as_adjacency_matrix(cur_sbm))

  # Replicate M time for current theta_ value
  sim_results <- lapply(seq(M), function(m){

    # Sample treatment vector (unit level vector but cluster randomization)
    cur_Z <- Z_ber_clusters(N_clusters = N_clusters,
                            N_each_cluster_vec = N_each_cluster_vec,
                            p = pz_cluster)
    # Compute exposures
    cur_expos <- generate_exposures_threshold(A = sbm_adj_mat,
                                              Z = cur_Z,
                                              threshold = threshold)

    # Compute observed outcomes using adj mat from SBM
    cur_outcomes <- generate_po(cur_expos, base.po.noise)


    # Estimate ATE using HT
    n = sum(N_each_cluster_vec)
    cur_ate = n^{-1}*2*(sum(cur_outcomes[cur_Z==1]) - sum(cur_outcomes[cur_Z==0]))


    expos_table = paste0(names(table(cur_expos)),"=",table(cur_expos), collapse = "; ")

    data.table(y_hat = cur_ate,
               y_true = 1, # Oracle ATE
               theta_ = theta_,
               expos_table = expos_table,
               m = m)
  })

  return(rbindlist(sim_results))
}


Bias_of_clusters_ATE <- function(Nv,
                                 N_clusters,
                                 N_each_cluster_vec,
                                 theta_seq,
                                 pz_cluster,
                                 threshold,
                                 M,
                                 distance_matrix=NULL,
                                 base.po.noise = NULL){


  full_results <- mclapply(theta_seq, function(theta_){

    # Estimate ATE for current theta_
    Cluster_randomization_CE_estimation(N_clusters = N_clusters,
                                        N_each_cluster_vec = N_each_cluster_vec,
                                        theta_ = theta_,
                                        Nv = Nv,
                                        pz_cluster = pz_cluster,
                                        threshold = threshold,
                                        M = M,
                                        distance_matrix = distance_matrix,
                                        base.po.noise = base.po.noise)

  },
  mc.cores = 1)
  # mc.cores = availableCores()/2)


  return(rbindlist(full_results))
}



create_noisy_adj_mat <- function(true_adj_mat,
                                 n_v,
                                 alpha_,
                                 beta_){

  # Save only lower triangular matrix

  lower_tri_idx <- lower.tri(true_adj_mat)

  noisy_adj_mat <- true_adj_mat

  noisy_adj_mat[!lower_tri_idx] <- NA

  # Remove existing ties with prob. beta_
  # Add ties with prob. alpha_

  ind_ties <- which(noisy_adj_mat==1)
  ind_no_ties <- which(noisy_adj_mat==0)
  noisy_adj_mat[ind_ties] <- rbinom(length(ind_ties),1,1-beta_)
  noisy_adj_mat[ind_no_ties] <- rbinom(length(ind_no_ties),1,alpha_)

  # Compute the full noisy matrix from the lower tri one

  # diag is 0
  diag(noisy_adj_mat) <- 0

  # Force the upper tri be the same as the lower tri

  noisy_adj_mat <- as.matrix(Matrix::forceSymmetric(noisy_adj_mat,uplo="L"))

  return(noisy_adj_mat)
}



create_censored_adj_mat <- function(true_adj_mat, K){
  censored.A <- true_adj_mat
  # Compute degrees
  deg.A <- rowSums(true_adj_mat)
  # Find which nodes with degree > K (and use permutation of it)
  to.censored.idx <- permute(which(deg.A > K))
  for (idx in to.censored.idx) {
    # Get current `idx` peers
    peers <- censored.A[,idx]
    n.peers <- sum(peers)
    if(n.peers > K){
      # sample n.peers-K units and remove the edge
      peers.ties.idx <- which(peers == 1)
      peers.to.remove <- sample(peers.ties.idx, n.peers - K, F)
      peers[peers.to.remove] <- 0
      # Update adj. matrix
      censored.A[,idx] <- peers
      censored.A[idx,] <- peers
    }
  }
  return(censored.A)
}



Bias_of_Noisy_Network_Simulation <- function(Nv,
                                             pz,
                                             p_v,
                                             net_model_name,
                                             # K_vec,
                                             beta_vec,
                                             threshold,
                                             adj_true_network,
                                             M,
                                             R,
                                             exposures_vec,
                                             exposures_contrast,
                                             base.po.noise = NULL,
                                             compute.exact.bias = TRUE){
  # Run the simulation
  full_results <- mclapply(beta_vec, function(beta){

    # Create noisy network out of censored adj mat
    if(beta>0){
      noisy_adj_mat <- create_noisy_adj_mat(true_adj_mat = adj_true_network,
                                            n_v =  Nv,
                                            alpha_ = beta/100,
                                            beta_ = beta)
    } else{
      noisy_adj_mat <- adj_true_network
    }

    jaccard_idx <- jaccard_edgeset_similarity(adj_true_network,noisy_adj_mat)

    One_network_CE_estimation(true_adj_mat = adj_true_network,
                              noisy_adj_mat = noisy_adj_mat,
                              network_model = net_model_name,
                              M = M,
                              Nv = Nv,
                              pz = pz,
                              R = R,
                              # param = K,
                              param = beta,
                              # alpha_ = alph_bet[1],
                              # beta_ = alph_bet[2],
                              threshold = threshold,
                              exposures_vec = exposures_vec,
                              exposures_contrast = exposures_contrast,
                              base.po.noise = base.po.noise,
                              jaccard_idx = jaccard_idx,
                              compute.exact.bias = compute.exact.bias)


  },
  mc.cores = 1)
  # mc.cores = detectCores())


  return(rbindlist(full_results))
}

Bias_of_Censored_Network_Simulation <- function(Nv,
                                                pz,
                                                p_v,
                                                net_model_name,
                                                K_vec,
                                                # beta_vec,
                                                threshold,
                                                adj_true_network,
                                                M,
                                                R,
                                                exposures_vec,
                                                exposures_contrast,
                                                base.po.noise = NULL,
                                                compute.exact.bias = TRUE){
  # Run the simulation
  full_results <- mclapply(K_vec, function(K){

    # Create censored network
    if(K > 0){
      censored_adj_mat <- create_censored_adj_mat(true_adj_mat = adj_true_network,
                                                  K = K)
    } else{
      censored_adj_mat <- adj_true_network
    }

    jaccard_idx <- jaccard_edgeset_similarity(adj_true_network,censored_adj_mat)

    One_network_CE_estimation(true_adj_mat = adj_true_network,
                              noisy_adj_mat = censored_adj_mat,
                              network_model = net_model_name,
                              M = M,
                              Nv = Nv,
                              pz = pz,
                              R = R,
                              param = K,
                              # param = beta,
                              # alpha_ = alph_bet[1],
                              # beta_ = alph_bet[2],
                              threshold = threshold,
                              exposures_vec = exposures_vec,
                              exposures_contrast = exposures_contrast,
                              base.po.noise = base.po.noise,
                              jaccard_idx = jaccard_idx,
                              compute.exact.bias = compute.exact.bias)


  },
  mc.cores = 1)
  # mc.cores = detectCores())


  return(rbindlist(full_results))
}



bias_variance_full_iter_in_setup <- function(R,
                                             k,
                                             iter,
                                             Pz_function,
                                             pz,
                                             true.adj.mat,
                                             A.list,
                                             adj.mat.names,
                                             with.adj,
                                             exposures_contrast, exposures_vec,
                                             threshold,
                                             base.po.noise=NULL){
  # For specific `k` value, A.list (subset of original), and other parameters,
  # Compute Prob. matrices and run `iter` iteration

  n <- nrow(A.list[[1]])

  # Need to compute probabilities matrices only once
  Prob_mat_list <- Get_prob_matrices_list(R = R, n = n,
                                          Pz_function = Pz_function,
                                          pz_func_args = list(n=n, p=pz),
                                          A.list = A.list,
                                          exposures_contrast = exposures_contrast,
                                          exposures_vec = exposures_vec,
                                          exposure_func = generate_exposures_threshold,
                                          # threshold = threshold)
                                          exposure_func_args = list(threshold))

  # `iter` times sample new Z, generate PO using true.adj.mat and estimate
  # using MR estimator and adj. matrices of A.list

  multiple_iter_results <- mclapply(seq(iter), function(m){
    # Generate data
    cur.Z <- Pz_function(n=n, p=pz)
    true.expos <- generate_exposures_threshold(A = true.adj.mat,
                                               Z = cur.Z,
                                               threshold = threshold)
    Y.obs <- generate_po(true.expos, base.po.noise)

    # Estimate CE using MR
    CE.dt <- rbindlist(MR_CE_estimator(Z.obs = cur.Z,
                                       Y.obs = Y.obs,
                                       A.list = A.list,
                                       exposures_contrast = exposures_contrast,
                                       exposures_vec = exposures_vec,
                                       Prob_matrices_list = Prob_mat_list,
                                       exposure_func = generate_exposures_threshold,
                                       exposure_func_args = list(threshold)))
    # Update results DT
    CE.dt[,`:=`(ce_contrast = sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])}),
                iter = m,
                K = k,
                with.true.adj = with.adj,
                adj.mat.used = paste0(adj.mat.names,collapse = ", ")
    )]

  },
  mc.cores = 1
  # mc.cores = availableCores()/2)
  # detectCores()/2
  )

  return(rbindlist(multiple_iter_results))

}

bias_variance_tradeoff_one_k_fullrun <- function(R,
                                                 k,
                                                 iter,
                                                 Pz_function, pz,
                                                 true.adj.mat,
                                                 true.adj.name,
                                                 A.list,
                                                 exposures_contrast, exposures_vec,
                                                 threshold,
                                                 base.po.noise = NULL){
  # Run bias-variance trade off simulation for all k-length adj. matrices combination

  names_combin <- combinations(n = length(A.list),
                               # r = k,
                               r = length(A.list),
                               v = names(A.list))

  one_k_results <- apply(names_combin, 1, function(cur.names){

    cur.a.lst <- A.list[cur.names]
    with.adj <- (true.adj.name %in% cur.names)

    bias_variance_full_iter_in_setup(R = R,
                                     k = k,
                                     iter = iter,
                                     Pz_function = Pz_function,
                                     pz = pz,
                                     true.adj.mat = true.adj.mat,
                                     A.list = cur.a.lst,
                                     adj.mat.names = cur.names,
                                     with.adj = with.adj,
                                     exposures_contrast = exposures_contrast,
                                     exposures_vec = exposures_vec,
                                     threshold = threshold,
                                     base.po.noise = base.po.noise)

  })

  return(rbindlist(one_k_results))
}

Bias_variance_tradeoff_in_MR <- function(R,
                                         K,
                                         iter,
                                         Pz_function, pz,
                                         true.adj.mat,
                                         true.adj.name,
                                         A.list,
                                         exposures_contrast, exposures_vec,
                                         threshold,
                                         base.po.noise = NULL){
  # Run bias-variance simulations for the MR estimator

  bias_variance_results <- lapply(seq(K), function(k){

    bias_variance_tradeoff_one_k_fullrun(R = R,
                                         k = k,
                                         iter = iter,
                                         Pz_function = Pz_function,
                                         pz = pz,
                                         true.adj.mat = true.adj.mat,
                                         true.adj.name = true.adj.name,
                                         A.list = A.list,
                                         exposures_contrast = exposures_contrast,
                                         exposures_vec = exposures_vec,
                                         threshold = threshold,
                                         base.po.noise = base.po.noise)
  })

  return(rbindlist(bias_variance_results))

}


jaccard_edgeset_similarity <- function(A1, A2) {
  # Transform to graph object
  G1 <- graph_from_adjacency_matrix(A1)
  G2 <- graph_from_adjacency_matrix(A2)
  # Compute intersection and union of edges
  inter <- length(E(G1 %s% G2))
  un <- length(E(G1 %u% G2))
  if (un == 0) {
    0
  } else {
    inter/un
  }
}

generate_Q_matrix <- function(N_clusters,
                              within_prob_vec,
                              between_prob_vec){
  # Generate Q matrix of SBM
  Q <- matrix(NA, N_clusters, N_clusters)
  diag(Q) <- within_prob_vec
  Q[lower.tri(Q)] <- between_prob_vec
  Q <- as.matrix(Matrix::forceSymmetric(Q, uplo = "L"))
  return(Q)
}

generate_clusters_contamination <- function(N_units,
                                             N_clusters,
                                             N_each_cluster_vec,
                                             within_prob_vec,
                                             between_prob_vec){
  # Generate network with cross-clusters contamination
  # Q.mat dimension are N_clusters*N_clusters
  Q.mat <- generate_Q_matrix(N_clusters = N_clusters,
                             within_prob_vec = within_prob_vec,
                             between_prob_vec = between_prob_vec)
  # SBM adj mat dimension are N_units*N_units
  SBM.adj.mat <- igraph::as_adjacency_matrix(igraph::sample_sbm(n = N_units,
                                                               pref.matrix = Q.mat,
                                                               block.sizes = N_each_cluster_vec,
                                                               directed = FALSE))
  return(SBM.adj.mat)
}

Single_PBA_CRT_iteration <- function(N_units,
                                     N_clusters,
                                     N_each_cluster_vec,
                                     prior_func,
                                     prior_func_args,
                                     between_prob_func,
                                     X.obs,
                                     Z.obs,
                                     Y.obs,
                                     Pz_function,
                                     pz_func_args,
                                     exposure_mapping = generate_exposures_threshold){
  # Sample theta (vector or singleton)
  n_theta <- length(prior_func)
  if(n_theta==1){
    cur.theta <- do.call(prior_func, prior_func_args)
  } else{
    cur.theta <- vector("numeric",n_theta)
    for (i in seq(n_theta)){
      cur.theta[i] <- do.call(prior_func[[i]], prior_func_args[[i]])
    }
  }
  # Generate cross-clusters edges probs.
  if(is.null(between_prob_func)){ # Same prob for all clusters
    between_prob_vec <- cur.theta
  } else{ # custom prob. using covariates X.obs 
    between_prob_vec <- do.call(between_prob_func, list(cur.theta, X.obs))
  }
  # Generate contaminated network
  cur.network <- generate_clusters_contamination(N_units = N_units,
                                                 N_clusters = N_clusters,
                                                 N_each_cluster_vec = N_each_cluster_vec,
                                                 within_prob_vec = rep(1,N_clusters),
                                                 between_prob_vec = between_prob_vec)
  
  # Get prob. matrices for current censored network
  prob.mat <- Get_prob_matrices_list(R = 10^4,
                                     n = N_units,
                                     Pz_function = Pz_function,
                                     pz_func_args = pz_func_args,
                                     A.list = list(cur.network),
                                     exposures_contrast = list(c("c11","c00")),
                                     exposures_vec = c("c11","c00"),
                                     threshold = rep(0,N_units))
  
  # Estimate causal effects
  # CE.estimate <- rbindlist(MR_CE_estimator(Z.obs = Z.obs,
  #                                    Y.obs = Y.obs,
  #                                    A.list = list(cur.network),
  #                                    exposures_contrast = list(c("c11","c00")),
  #                                    exposures_vec = c("c11","c00"),
  #                                    Prob_matrices_list = prob.mat,
  #                                    threshold = rep(0,N_units)))
  # 
  CE.estimate <- MR_CE_estimator(Z.obs = Z.obs,
                                     Y.obs = Y.obs,
                                     A.list = list(cur.network),
                                     exposures_contrast = list(c("c11","c00")),
                                     exposures_vec = c("c11","c00"),
                                     Prob_matrices_list = prob.mat,
                                     threshold = rep(0,N_units))
  
  
  CE.estimate <- data.table(do.call(rbind,CE.estimate))
  
  
  # Update results DT
  # CE.estimate <- NMR_estimator(A.list = list(A=cur.network),
  #                              Z.obs = Z.obs,
  #                              Y.obs = Y.obs,
  #                              Pz_function = Pz_function,
  #                              pz_func_args = pz_func_args,
  #                              exposures_vec = c("c11","c00"),
  #                              exposures_contrast = list(c("c11","c00")),
  #                              exposure_func = generate_exposures_threshold,
  #                              exposure_func_args = list(threshold = rep(0,N_units)))
  
  # Add random error via normal approximation
  CE.estimate$ht_ce_w_re <- CE.estimate$ht_ce + rnorm(1,0,sqrt(CE.estimate$var_ht_ce))
  CE.estimate$hajek_ce_w_re <- CE.estimate$hajek_ce + rnorm(1,0,sqrt(CE.estimate$var_hajek_ce))
  
  CE.estimate$ce_contrast <- sapply(list(c("c11","c00")), function(x){paste0(x[1],"-",x[2])})
  
  return(CE.estimate)
}

PBA_for_CRT <- function(N_units,
                        N_clusters,
                        N_each_cluster_vec,
                        N_iterations,
                        prior_func,
                        prior_func_args,
                        between_prob_func = NULL,
                        X.obs = NULL,
                        Z.obs,
                        Y.obs,
                        Pz_function,
                        pz_func_args){
  
  PBA.results <- mclapply(seq(N_iterations),function(m){
        cur.results <- Single_PBA_CRT_iteration(N_units = N_units,
                                                N_clusters = N_clusters,
                                                N_each_cluster_vec = N_each_cluster_vec,
                                                prior_func = prior_func,
                                                prior_func_args = prior_func_args,
                                                between_prob_func = between_prob_func,
                                                X.obs = X.obs,
                                                Z.obs = Z.obs, 
                                                Y.obs = Y.obs,
                                                Pz_function = Pz_function,
                                                pz_func_args = pz_func_args)
        cur.results$iter <- m
        cur.results
  },
  mc.cores = 1
  # mc.cores = 8
  # mc.cores = parallel::detectCores()/2
  )
  return(rbindlist(PBA.results))
}

generate_perturbed_network <- function(N_units,
                                       A.sp,
                                       edges_prob_func,
                                       prior_func_list,
                                       prior_func_args_list,
                                       X.obs = NULL){
  # Sample from prior
  theta.length <- length(prior_func_list)
  theta.vec <- vector("numeric",theta.length)
  for (i in seq(theta.length)){
    theta.vec[i] <- do.call(prior_func_list[[i]],prior_func_args_list[[i]])
  }
  # Generate P_\theta prob matrix (lower tri matrix with upper tri = 0 or NA)
  if(is.null(X.obs)){
    prob.lower.tri.matrix <- do.call(edges_prob_func, list(theta.vec, A.sp))
  } else{
    prob.lower.tri.matrix <- do.call(edges_prob_func, list(theta.vec, A.sp, X.obs))
  }
  # Sample edges
  edges.prob.vec <- prob.lower.tri.matrix[lower.tri(prob.lower.tri.matrix)]
  edges.lower.tri <- rbinom(n = length(edges.prob.vec), size = 1, prob = edges.prob.vec)
  # Save as A.star
  A.star <- matrix(NA, N_units, N_units)
  A.star[lower.tri(A.star)] <- edges.lower.tri
  diag(A.star) <- 0
  # Make it symmetric (undirected)
  A.star <- Matrix::forceSymmetric(A.star, uplo = "L")
  return(A.star)
}

Single_PBA_general_iteration <- function(N_units,
                                         A.sp,
                                         edges_prob_func,
                                         prior_func_list,
                                         prior_func_args_list,
                                         Z.obs,
                                         Y.obs,
                                         X.obs = NULL,
                                         Pz_function,
                                         pz_func_args,
                                         exposures_vec,
                                         exposures_contrast,
                                         exposures_thresholds = NULL){
    
  # Generate A.star (perturbed version of A.sp)
  cur.A.star <- generate_perturbed_network(N_units = N_units,
                                           A.sp = A.sp,
                                           edges_prob_func = edges_prob_func,
                                           prior_func_list = prior_func_list,
                                           prior_func_args_list = prior_func_args_list,
                                           X.obs = X.obs)
  
  
  if(is.null(exposures_thresholds)){exposures_thresholds = rep(0, N_units)}
  
  # Get prob. matrices for current censored network
  prob.mat <- Get_prob_matrices_list(
                                      # R = 10^4,
                                      R = 10,
                                     n = N_units,
                                     Pz_function = Pz_function,
                                     pz_func_args = pz_func_args,
                                     A.list = list(cur.A.star),
                                     exposures_contrast = exposures_contrast,
                                     exposures_vec = exposures_vec,
                                     threshold = exposures_thresholds)
  

  cur.CE.estimates <- MR_CE_estimator(Z.obs = Z.obs,
                                           Y.obs = Y.obs,
                                           A.list = list(cur.A.star),
                                           exposures_contrast = exposures_contrast,
                                           exposures_vec = exposures_vec,
                                           Prob_matrices_list = prob.mat,
                                           threshold = exposures_thresholds)
  
  cur.CE.estimates <- data.table(do.call(rbind,cur.CE.estimates))
  
  # Add random error via normal approximation
  cur.CE.estimates$ht_ce_w_re <- mapply(rnorm,
                                        n=1,
                                        mean=cur.CE.estimates$ht_ce,
                                        sd=sqrt(cur.CE.estimates$var_ht_ce))
  
  cur.CE.estimates$hajek_ce_w_re <- mapply(rnorm,
                                        n=1,
                                        mean=cur.CE.estimates$hajek_ce,
                                        sd=sqrt(cur.CE.estimates$var_hajek_ce))
  
  cur.CE.estimates$ce_contrast <- sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])})
  
  return(cur.CE.estimates)
}

PBA_general <- function(N_units,
                        N_iterations,
                        A.sp,
                        edges_prob_func,
                        prior_func_list,
                        prior_func_args_list,
                        Z.obs,
                        Y.obs,
                        X.obs = NULL,
                        Pz_function,
                        pz_func_args,
                        exposures_vec,
                        exposures_contrast,
                        exposures_thresholds = NULL){
  
  PBA.results <- mclapply(seq(N_iterations),function(m){
    CE.estimates <- Single_PBA_general_iteration(N_units = N_units,
                                                 A.sp = A.sp,
                                                 edges_prob_func = edges_prob_func,
                                                 prior_func_list = prior_func_list,
                                                 prior_func_args_list = prior_func_args_list,
                                                 Z.obs = Z.obs,
                                                 Y.obs = Y.obs,
                                                 X.obs = X.obs,
                                                 Pz_function = Pz_function,
                                                 pz_func_args = pz_func_args,
                                                 exposures_vec = exposures_vec,
                                                 exposures_contrast = exposures_contrast,
                                                 exposures_thresholds = exposures_thresholds)
    CE.estimates[,iter := m]
    CE.estimates
  },
  mc.cores = 1
  # mc.cores = 8
  # mc.cores = parallel::detectCores()/2
  )
  return(rbindlist(PBA.results))
}


# Data analysis -----------------------------------------------------------

compute_prob_indicator_palluck <- function(R, n, Pz_function, pz_func_args, A.list, exposures_vec){
  # A function that returns a list of nXr indicator matrices (one for each exposures),
  # that indicate if unit i in sample r is exposed to some ck under all A \in A.list
  # The output is used for computation of exposures probabilities


  # Init a list of nXR indicator matrices
  Ind.mat.list <- vector("list",length(exposures_vec))
  for (i in seq(length(exposures_vec))) {
    assign(paste0("I.",exposures_vec[i]), matrix(NA, nrow = n, ncol = R))
    Ind.mat.list[[i]] <- get(paste0("I.",exposures_vec[i]))
  }
  names(Ind.mat.list) <- exposures_vec

  # Sample Z ~ PZ R times and update the indicators matrices
  for (r in seq(R)) {
    # cur.Z <- Pz_function(n=n, p=pz, schid=schid)
    cur.Z <- do.call(Pz_function, pz_func_args)
    cur.exposure.mat <- exposures_under_network_list_palluck(Z=cur.Z$Z.vec,
                                                             school.treated = cur.Z$school.treated,
                                                             A.list=A.list)
    for (ck in exposures_vec) {
      # Update indiactor matrices for each exposures and r
      Ind.mat.list[[ck]][,r] <- specific_exposure_under_network_list(exposures_mat=cur.exposure.mat,
                                                                     ck=ck)
    }
  }
  return(Ind.mat.list)
}


exposures_under_network_list_palluck <- function(Z, school.treated, A.list){
  # for m adj. matrices in A.list and n units, this function compute the nXm matrix of exposures under
  # each adj. matrix for treatment vector Z
  return(
    sapply(A.list, function(A){compute.exposures.palluck(Z, school.treated, A)})
  )}


treat.randomization.palluck <- function(n, pz=0.5, schid){
  trt.vector <- vector("numeric", length(schid))
  uniq.schools <- unique(schid)
  # get treated schools
  trt.schools <- sample(x = uniq.schools,size = round(length(uniq.schools)/2), replace = FALSE)
  # sample unit level treatment for each treated school
  for (schl in trt.schools) {
    trt.vector[schid == schl] <- rbinom(n=sum(schid==schl), size = 1, prob = pz)
  }
  is.treated.school <- schid %in% trt.schools
  return(list(Z.vec = trt.vector, school.treated = is.treated.school))
}


compute.exposures.palluck <- function(Z, school.treated, adj.mat){
  # Compute num neighbors
  # Compute signs
  indirect_sign <- as.numeric((adj.mat %*% Z) > 0)
  deg <- rowSums(adj.matrix.either)
  # Save exposures
  exposures <- vector(length = nrow(adj.mat))

  exposures[as.logical(Z*indirect_sign*school.treated)] <- "c111"
  exposures[as.logical(Z*(1-indirect_sign)*school.treated)] <- "c101"
  exposures[as.logical((1-Z)*indirect_sign*school.treated)] <- "c011"
  exposures[as.logical((1-Z)*(1-indirect_sign)*school.treated)] <- "c001"
  exposures[as.logical(1-school.treated)] <- "c000"
  exposures[deg==0] <- "ISOLATED"

  return(exposures)
}

MR_CE_estimator_palluck <- function(Z.obs,school.treated ,Y.obs, A.list,
                                    exposures_contrast, exposures_vec,
                                    Prob_matrices_list,
                                    threshold,
                                    estimate.n = FALSE){

  expos.obs.mat <- exposures_under_network_list_palluck(Z = Z.obs,school.treated = school.treated,
                                                        A.list = A.list)

  ce.results <- vector("list", length(exposures_contrast))
  names(ce.results) <- sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])})

  for (expos in exposures_contrast) {
    ck <- expos[1]
    cl <- expos[2]
    cntrast_name <- paste0(ck,"-",cl)
    P.k <- Prob_matrices_list[[ck]]
    P.l <- Prob_matrices_list[[cl]]
    P.kl <- Prob_matrices_list[[cntrast_name]]
    ce_estimate <- MR_CE_wrapper(Y.obs = Y.obs,
                                 expos.obs.mat = expos.obs.mat,
                                 ck = ck,
                                 cl = cl,
                                 P.k = P.k,
                                 P.l = P.l,
                                 P.kl = P.kl,
                                 estimate.n = estimate.n)
    ce.results[[cntrast_name]] <- ce_estimate
  }
  return(ce.results)
}




