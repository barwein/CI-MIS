###
# Script that analyze the CRT in Kenya conducted by Venturo-Conerly et al., 2022
###


# Load data and libraries -------------------------------------------------


library(data.table)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(ggtext)
library(kableExtra)


source("Reproducibility/Estimation_functions.R")

library(MisspecifiedInterference)

school.df <- fread("Reproducibility/Data_analyses/Kenya_CRT/imputed_data.csv")

# Clean Data --------------------------------------------------------------


# length(unique(school.df$Participant_ID)) # 895 particpants from 2 schools and
# length(unique(school.df$Class)) # 24 classes (12 per school)
#
# table(school.df$Condition) # Study Skills = Control group
# table(school.df$Condition, school.df$School) # Study Skills = Control group

school.df.pre <- school.df[Time==0,]
school.df.post <- school.df[Time==1,]

setkeyv(school.df.pre, "Participant_ID")
setkeyv(school.df.post, "Participant_ID")

# Combine pre- and post- data
school.relevant.only <- school.df.pre[school.df.post]

school.relevant.only <- school.relevant.only[,.(Class, School, Participant_ID,
                                                Gender,
                                                Condition,
                                                GAD_Total, i.GAD_Total,
                                                PHQ_Total, i.PHQ_Total)]

# Compute difference in outcome (post-pre)
school.relevant.only[,`:=`(gad.diff = i.GAD_Total - GAD_Total,
                           phq.diff = i.PHQ_Total - PHQ_Total)]

# Define new version of treatment
school.relevant.only[,Treatment := as.numeric(Condition %in% c("Values","Growth"))]

# Save DF

write.csv(school.relevant.only,
          "Reproducibility/Data_analyses/Kenya_CRT/school.df.csv",
          row.names = F)

# ATE under the null ------------------------------------------------------

# Estimate ATE using HT estimator (p(T=1)=0.5 for each unit thus weighting by 1/0.5 = 2)

n <- nrow(school.relevant.only)

gad.class.sum <- tapply(school.relevant.only$gad.diff, school.relevant.only$Class, sum)
gad.trt.classs <- tapply(school.relevant.only$Treatment, school.relevant.only$Class, function(x){all(x==1)})


# ATE and SE based on Aronow&Midelton 2013; they expressed cluster estimators using HT.

null.ATE.GAD <- n^{-1}*2*(sum(gad.class.sum[gad.trt.classs]) - sum(gad.class.sum[!gad.trt.classs]))

null.Var.ATE.GAD <- n^{-2}*(sum(gad.class.sum[gad.trt.classs]^2) - sum(gad.class.sum[!gad.trt.classs]^2))

null.gad.ci.low <- null.ATE.GAD - (1.96 * sqrt(null.Var.ATE.GAD))
null.gad.ci.high <- null.ATE.GAD + (1.96 * sqrt(null.Var.ATE.GAD))

paste0("GAD effect estimate [CI] is ", round(null.ATE.GAD,3),
       " [", round(null.gad.ci.low,3), ", ",
       round(null.gad.ci.high,3),"]")


# PHQ
phq.class.sum <- tapply(school.relevant.only$phq.diff, school.relevant.only$Class, sum)
phq.trt.classs <- tapply(school.relevant.only$Treatment, school.relevant.only$Class, function(x){all(x==1)})

null.ATE.phq <- n^{-1}*2*(sum(phq.class.sum[phq.trt.classs]) - sum(phq.class.sum[!phq.trt.classs]))

null.Var.ATE.phq <- n^{-2}*(sum(phq.class.sum[phq.trt.classs]^2) - sum(phq.class.sum[!phq.trt.classs]^2))

null.phq.ci.low <- null.ATE.phq - (1.96 * sqrt(null.Var.ATE.phq))
null.phq.ci.high <- null.ATE.phq + (1.96 * sqrt(null.Var.ATE.phq))

paste0("PHQ effect estimate [CI] is ", round(null.ATE.phq,3),
       " [", round(null.phq.ci.low,3), ", ",
       round(null.phq.ci.high,3),"]")


# PBA ----------------------------------------------------------------------

# Read data
kenya_df <- fread("Reproducibility/Data_analyses/Kenya_CRT/school.df.csv")

# Order by Class, School
setorderv(kenya_df, c("Class","School"))

school.names <- unique(kenya_df$School) # "KDSS" "MGSS"

n.each.KDSS <- c(table(kenya_df[School == "KDSS", Class]))
n.each.MGSS <- c(table(kenya_df[School == "MGSS", Class]))

N_clusters <- length(unique(kenya_df$Class))
N_each_cluster_vec <- c(n.each.KDSS, n.each.MGSS)

# Baseline GAD estimate tau(c11,c00)

baseline.cluster.net <- generate_clusters_contamination(N_units = n,
                                                        N_clusters = N_clusters,
                                                        N_each_cluster_vec = N_each_cluster_vec,
                                                        within_prob_vec = 1,
                                                        between_prob_vec = 0)
set.seed(151502)
baseline.GAD.estimate <- NMR_estimator(A.list = list(A=baseline.cluster.net),
                                       Z.obs = kenya_df$Treatment,
                                       Y.obs = kenya_df$gad.diff,
                                       Pz_function = Z_ber_clusters,
                                       pz_func_args = list(N_clusters = N_clusters,
                                                           N_each_cluster_vec = N_each_cluster_vec,
                                                           p = 0.5),
                                       exposures_vec = c("c11","c00"),
                                       exposures_contrast = list(c("c11","c00")),
                                       exposure_func = generate_exposures_threshold,
                                       exposure_func_args = list(threshold = rep(0,n)))

baseline.hajek <- round(baseline.GAD.estimate$hajek_ce,3)
# CI via normal approximation and conservative variance estimator
baseline.se.times.quantile <- round(sqrt(baseline.GAD.estimate$var_hajek_ce)*1.96,3)
baseline.point.ci <- paste0(baseline.hajek," [",
                            baseline.hajek - baseline.se.times.quantile,
                            ", ", baseline.hajek + baseline.se.times.quantile,"]")

# PBA aux functions

Between_cluster_edges_prob <- function(theta.vec, school.id){
  # Function that generates cross-clusters edges with prob theta_{ij} for units in the same school
  # school.id is cluster level covarite
  q <- length(school.id) # number of clusters
  theta.mat <- matrix(NA,q,q)
  theta.mat[lower.tri(theta.mat)] <- theta.vec # matrix of theta values
  between.prob.mat <- matrix(NA,q,q) # results matrix
  for (i in seq(q-1)){
    for (j in seq(i+1,q)) {
        # P(A_{ij}=1) = I(school_i=school_j)*theta_ij 
        # if theta.vec is singelton then theta_ij = theta 
        between.prob.mat[j,i] <- theta.mat[j,i]*(school.id[i]==school.id[j])    
    }
  }
  return(between.prob.mat[lower.tri(between.prob.mat)])
}


Edges_prob_by_school_and_gender <- function(theta.vec, A.sp, X){
  # Function that generates cross-clusters edges with prob that depends on both school and gender
  school <- X[[1]]
  gender <- X[[2]]
  N <- dim(A.sp)[1]
  prob.mat <- matrix(NA, N, N) 
  diag(prob.mat) <- 0
  for (i in seq(N-1)) {
    for (j in seq(i+1,N)) {
      # P(A_ij=1 | A^sp_ij = 0) = I(school_i=school_j){theta[1] + theta[2]*(gender_i=gender_j)}
      # P(A_ij=1 | A^sp_ij = 1) = 1
      prob.mat[j,i] <- ifelse(A.sp[j,i] == 1,
                              1,
                              (school[i]==school[j])*(theta.vec[1] + theta.vec[2]*(gender[i]==gender[j])))
    }
  }
  return(prob.mat)
}

# Same theta
# Uniform(0,0.005) prior
set.seed(62619)
PBA.same.theta.uniform.prior <- PBA_for_CRT(N_units = n,
                                N_clusters = N_clusters,
                                N_each_cluster_vec = N_each_cluster_vec,
                                N_iterations = 1e3,
                                prior_func = runif,
                                prior_func_args = list(n=1,min=0,max=0.005),
                                between_prob_func = Between_cluster_edges_prob,
                                X.obs = rep(c("KDSS","MGSS"),each=12),
                                Z.obs = kenya_df$Treatment,
                                Y.obs = kenya_df$gad.diff,
                                Pz_function = Z_ber_clusters,
                                pz_func_args = list(N_clusters = N_clusters,
                                                    N_each_cluster_vec = N_each_cluster_vec,
                                                    p = 0.5))
# # Beta(0.25,20) prior
set.seed(62620)
PBA.same.theta.beta.prior <- PBA_for_CRT(N_units = n,
                                N_clusters = N_clusters,
                                N_each_cluster_vec = N_each_cluster_vec,
                                N_iterations = 1e3,
                                prior_func = rbeta,
                                prior_func_args = list(n=1,shape1=0.25,shape2=25),
                                between_prob_func = Between_cluster_edges_prob,
                                X.obs = rep(c("KDSS","MGSS"),each=12),
                                Z.obs = kenya_df$Treatment,
                                Y.obs = kenya_df$gad.diff,
                                Pz_function = Z_ber_clusters,
                                pz_func_args = list(N_clusters = N_clusters,
                                                    N_each_cluster_vec = N_each_cluster_vec,
                                                    p = 0.5))

prior_func_multi <- rep(list(runif), N_clusters*(N_clusters-1)/2)
prior_func_args_multi <- rep(list(list(n=1,min=0,max=0.005)), N_clusters*(N_clusters-1)/2)

PBA.different.theta.uniform.prior <- PBA_for_CRT(N_units = n,
                                        N_clusters = N_clusters,
                                        N_each_cluster_vec = N_each_cluster_vec,
                                        N_iterations = 1e3,
                                        prior_func = prior_func_multi,
                                        prior_func_args = prior_func_args_multi,
                                        between_prob_func = Between_cluster_edges_prob,
                                        X.obs = rep(c("KDSS","MGSS"),each=12),
                                        Z.obs = kenya_df$Treatment,
                                        Y.obs = kenya_df$gad.diff,
                                        Pz_function = Z_ber_clusters,
                                        pz_func_args = list(N_clusters = N_clusters,
                                                            N_each_cluster_vec = N_each_cluster_vec,
                                                            p = 0.5))

# Different theta by school

# Uniform prior
prior_func_multi <- rep(list(runif), N_clusters*(N_clusters-1)/2)
prior_func_args_multi <- rep(list(list(n=1,min=0,max=0.005)), N_clusters*(N_clusters-1)/2)

set.seed(65130)
PBA.different.theta.uniform.prior <- PBA_for_CRT(N_units = n,
                                        N_clusters = N_clusters,
                                        N_each_cluster_vec = N_each_cluster_vec,
                                        N_iterations = 1e3,
                                        prior_func = prior_func_multi,
                                        prior_func_args = prior_func_args_multi,
                                        between_prob_func = Between_cluster_edges_prob,
                                        X.obs = rep(c("KDSS","MGSS"),each=12),
                                        Z.obs = kenya_df$Treatment,
                                        Y.obs = kenya_df$gad.diff,
                                        Pz_function = Z_ber_clusters,
                                        pz_func_args = list(N_clusters = N_clusters,
                                                            N_each_cluster_vec = N_each_cluster_vec,
                                                            p = 0.5))

# Beta prior
prior_func_multi <- rep(list(rbeta), N_clusters*(N_clusters-1)/2)
prior_func_args_multi <- rep(list(list(n=1,shape1=0.25,shape2=25)), N_clusters*(N_clusters-1)/2)

set.seed(65131)
PBA.different.theta.beta.prior <- PBA_for_CRT(N_units = n,
                                        N_clusters = N_clusters,
                                        N_each_cluster_vec = N_each_cluster_vec,
                                        N_iterations = 1e3,
                                        prior_func = prior_func_multi,
                                        prior_func_args = prior_func_args_multi,
                                        between_prob_func = Between_cluster_edges_prob,
                                        X.obs = rep(c("KDSS","MGSS"),each=12),
                                        Z.obs = kenya_df$Treatment,
                                        Y.obs = kenya_df$gad.diff,
                                        Pz_function = Z_ber_clusters,
                                        pz_func_args = list(N_clusters = N_clusters,
                                                            N_each_cluster_vec = N_each_cluster_vec,
                                                            p = 0.5))


# Base line Q is diag network (no cross clusters contamination)
Q.sp <- diag(1, N_clusters, N_clusters)
A.sp <- igraph::as_adjacency_matrix(igraph::sample_sbm(n = n,
                                                       pref.matrix = Q.sp,
                                                       block.sizes = N_each_cluster_vec,
                                                       directed = FALSE))
set.seed(651995)
PBA.theta.by.school.gender.uniform.prior <- PBA_general(N_units = n,
                                                      N_iterations = 2,
                                                      A.sp = A.sp,
                                                      edges_prob_func = Edges_prob_by_school_and_gender,
                                                      prior_func_list = rep(list(runif),2),
                                                      prior_func_args_list = list(list(n=1,min=0,max=0.003),
                                                                                  list(n=1,min=0,max=0.002)),
                                                      Z.obs = kenya_df$Treatment,
                                                      Y.obs = kenya_df$gad.diff,
                                                      X.obs = list(school = kenya_df$School,
                                                                   gender = kenya_df$Gender),
                                                      Pz_function = Z_ber_clusters,
                                                      pz_func_args = list(N_clusters = N_clusters,
                                                                          N_each_cluster_vec = N_each_cluster_vec,
                                                                          p = 0.5),
                                                      exposures_vec = c("c11","c00"),
                                                      exposures_contrast = list(c("c11","c00")))

# Point mass theta distribution

theta_grid <- seq(0.0005,0.005,0.0005)

set.seed(62619)
for (i in seq(length(theta_grid))){
  Kenya.PBA <- PBA_for_CRT(N_units = n,
                           N_clusters = N_clusters,
                           N_each_cluster_vec = N_each_cluster_vec,
                           N_iterations = 200,
                           prior_func = sample,
                           prior_func_args = list(x=theta_grid[i],p=1),
                           between_prob_func = Between_cluster_edges_prob,
                           X.obs = rep(c("KDSS","MGSS"),each=12),
                           Z.obs = kenya_df$Treatment,
                           Y.obs = kenya_df$gad.diff, 
                           Pz_function = Z_ber_clusters,
                           pz_func_args = list(N_clusters = N_clusters,
                                               N_each_cluster_vec = N_each_cluster_vec,
                                               p = 0.5))
  write.csv(x = Kenya.PBA,
            file = paste0("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/pointmass/PBA_results/Kenya_PBA_same_theta_point_mass_theta", i, ".csv"),
            row.names = FALSE)
}





# PBA graphics ------------------------------------------------------------

# Read files
Same.theta.uniform.prior <- fread("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_PBA_same_theta_uniform_prior.csv")
Same.theta.uniform.prior[, scenario := "same.theta.uniform"]

Same.theta.beta.prior <- fread("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_PBA_same_theta_beta_prior.csv")
Same.theta.beta.prior[, scenario := "same.theta.beta"]

different.theta.uniform.prior <- fread("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_PBA_different_theta_uniform_prior.csv")
different.theta.uniform.prior[, scenario := "different.theta.uniform"]

different.theta.beta.prior <- fread("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_PBA_different_theta_beta_prior.csv")
different.theta.beta.prior[, scenario := "different.theta.beta"]

By.school.gender.uniform <- fread("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_PBA_theta_by_school_gender_uniform_prior.csv")
By.school.gender.uniform[, scenario := "by.school.gender"]

# combine all
PBA.combined <- rbindlist(list(Same.theta.uniform.prior,
                               Same.theta.beta.prior,
                               different.theta.uniform.prior,
                               different.theta.beta.prior, 
                               By.school.gender.uniform))


PBA.combined <- melt.data.table(PBA.combined,
                                id.vars = c("iter","scenario"),
                                measure.vars = c("ht_ce","hajek_ce","ht_ce_w_re","hajek_ce_w_re"),
                                variable.name = "estimator")

ggplot(PBA.combined[estimator == "hajek_ce",], aes(x=value)) + 
  geom_histogram() +
  facet_wrap(~scenario, nrow = 1, scales = "free")

ggplot(PBA.combined[estimator == "hajek_ce_w_re",], aes(x=value)) + 
  geom_histogram() +
  facet_wrap(~scenario, nrow = 1, scales = "free")


# compute summary
PBA.summarized <- PBA.combined[,.(mean.esti = mean(value, na.rm=T),
                                  median.esti = median(value, na.rm=T),
                                  q025.esti = quantile(value, 0.025, na.rm=T),
                                  q975.esti = quantile(value, 0.975, na.rm=T)),
                               by = c("scenario","estimator")]

PBA.summarized[, with_re := ifelse(grepl("_w_re", PBA.summarized$estimator, fixed=T),"Yes","No")]
PBA.summarized[, estimator.type := ifelse(grepl("ht", PBA.summarized$estimator, fixed=T),"HT","Hajek")]

# add null results
PBA.summarized <- rbindlist(list(PBA.summarized,
                                 data.table(scenario = "null",
                                            estimator = "null",
                                            mean.esti = null.ATE.GAD,
                                            # mean.esti = baseline.hajek,
                                            median.esti = NA,
                                            # q025.esti = baseline.hajek - baseline.se.times.quantile,
                                            q025.esti = null.gad.ci.low,
                                            # q975.esti = baseline.hajek + baseline.se.times.quantile,
                                            q975.esti = null.gad.ci.high,
                                            with_re = "null",
                                            estimator.type = "null")))

PBA.summarized[,interval.length := q975.esti - q025.esti]

# Graphics

PBA.summarized$scenario <- factor(PBA.summarized$scenario, 
                                   levels = c("null",
                                              "same.theta.uniform",
                                              "same.theta.beta",
                                              "different.theta.uniform",
                                              "different.theta.beta",
                                              "by.school.gender"),
                                   labels = c("Baseline",
                                              "(I) & Uniform",
                                              "(I) & Beta",
                                              "(II) & Uniform",
                                              "(II) & Beta",
                                              "(III)"))

PBA.summarized$width <- ifelse(PBA.summarized$scenario == "Baseline", 0.2,0.4)

PBA.CI.figure <- ggplot(PBA.summarized[estimator.type %in% c("Hajek","null"),],
                        aes(x=scenario, color = with_re, group = with_re)) +
                  geom_errorbar(aes(ymin=q025.esti,
                                    ymax=q975.esti,
                                    width = width),
                                linewidth = 3,
                                alpha = 0.9,
                                position = position_dodge(0.4),
                                # width = 0.2,
                                show.legend = F) +
                  geom_point(aes(y=mean.esti),
                             size = 12,
                             alpha = .95,
                             position = position_dodge(0.4),
                             show.legend = F) +
                  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.35, linewidth = .8) +
                  scale_y_continuous(breaks = seq(-5,4,1)) +
                  scale_color_manual(values = c("null" = "#1D8A99","No" = "grey45","Yes"="grey25")) +
                  labs(x = "", y = "") + 
                  # theme_pubclean() +
                  theme_pubr() +
                  coord_flip() +
                  theme(axis.text.x = element_text(size =32, face = "bold"),
                        axis.text.y = element_text(size =32, face = "bold"), 
                        axis.ticks.y = element_blank())

ggsave(filename = "Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_CRT_PBA_CI_figure.jpeg",
       plot = PBA.CI.figure,
       width = 20,
       height = 12)

PBA.summarized[,point.and.CI := paste0(round(mean.esti,3),
                                      " [",round(q025.esti,3),", ",
                                      round(q975.esti,3),"]")]

PBA.summarized[with_re == "null", with_re := "Yes"]

PBA.summarized.casted <- dcast.data.table(PBA.summarized[estimator.type %in% c("Hajek","null")],
                          formula = "scenario ~ with_re", value.var = "point.and.CI")

kable(PBA.summarized.casted, format = "latex", booktabs=T,
              # col.names = c("Scenario","Without random error","With random error")) 
              col.names = c("Scenario","No","Yes")) %>%
              add_header_above(c(" " = 1, "Random error" = 2))

###
# Point mass \theta grid plots
###

# Merge to one CSV file

# csv_names <- dir("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/pointmass/")
# csv_names <- csv_names[c(1,seq(3,10),2)]
# path_name <- "Reproducibility/Data_analyses/Kenya_CRT/PBA_results/pointmass/"
# 
# theta_grid <- seq(0.0005,0.005,0.0005)
# theta_grid_perc <- paste0(round(theta_grid*100,2),"%")
# 
# combined_dt <- data.table()
# 
# for (i in seq(length(csv_names))) {
# 
#   curr_dt <- fread(paste0(path_name,csv_names[i]))
#   curr_dt[,param := theta_grid_perc[i]]
# 
#   combined_dt <- rbindlist(list(combined_dt,curr_dt))
# }
# 
# 
# write.csv(combined_dt,
#           "Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_PBA_same_theta_point_mass.csv",
#           row.names = FALSE)
# 

# Plot

PBA.point.mass.theta <- fread("Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_PBA_same_theta_point_mass.csv")

PBA.point.mass.theta.melt <- melt.data.table(PBA.point.mass.theta,
                                             id.vars = c("iter","param"),
                                             measure.vars = c("hajek_ce","hajek_ce_w_re"))


PBA.point.mass.boxplot <- ggplot(PBA.point.mass.theta.melt[variable=="hajek_ce",],
                                     aes(x=param, y=value)) +
                                geom_boxplot(fill = "grey45", alpha=0.7) +
                                geom_hline(yintercept = -0.658,
                                           color = "#1D8A99",
                                           linewidth = 1) + # Baseline point estimate 
                                labs(x=TeX("$\\theta$"),y="") +
                                scale_y_continuous(breaks = seq(-2.5,1,0.5)) +
                                theme_pubclean() +
                                theme(axis.text.x = element_text(size =28, face = "bold"),
                                      axis.text.y = element_text(size =26, face = "bold"),
                                      axis.title.x = element_text(size=38))

ggsave(filename = "Reproducibility/Data_analyses/Kenya_CRT/PBA_results/Kenya_CRT_PBA_pointmass_boxplot.jpeg",
       plot = PBA.point.mass.boxplot,
       width = 20,
       height = 12)






