###
# Script that analyze the social network field experiment (Palluck et al. 2016)
###

# Load libraries and files ------------------------------------------------

# Clean data, get networks, and define relevant variables

library(igraph)
library(Matrix)
library(stringr)
library(kableExtra)

source("Reproducibility/Estimation_functions.R")

library(MisspecifiedInterference)

options(na.action = na.omit)
# Read data
load("Reproducibility/Data_analyses/Palluck_et_al/37070-0001-Data.rda")
all.schools <- da37070.0001


# Clean data  -------------------------------------------------------------


## Disregard all observations for which both the school-id and the block-number of the school are unknown (1.14% of the observations, or 280 observations)
all.schools <- all.schools[!is.na(all.schools$SCHRB) & all.schools$ID != 999, ]

# Treatment variables to numeric
treat.to.numeric <- str_extract_all(all.schools$TREAT,"\\(?[0-9,.]+\\)?")
treat.to.numeric <- as.numeric(gsub("\\((.+)\\)","\\1",treat.to.numeric))
all.schools$TREAT.NUMERIC <- treat.to.numeric

schtreat.to.numeric <- str_extract(all.schools$SCHTREAT,"\\(?[0-9,.]+\\)?")
schtreat.to.numeric <- as.numeric(gsub("\\((.+)\\)","\\1",schtreat.to.numeric))
all.schools$SCHTREAT.NUMERIC <- schtreat.to.numeric

# indicator of eligible units
all.schools$ELIGIBLE <- 1*(all.schools$TREAT.NUMERIC>0)
all.schools$ELIGIBLE[is.na(all.schools$ELIGIBLE)] <- 0


# # Numeric version of OUTOFBLOCK
outofblock.to.numeric <- str_extract(all.schools$OUTOFBLOCK,"\\(?[0-9,.]+\\)?")
outofblock.to.numeric <- as.numeric(gsub("\\((.+)\\)","\\1",outofblock.to.numeric))
all.schools$OUTOFBLOCK.NUMERIC <- outofblock.to.numeric


# Remove units with OUTOFBLOCK=1 & BLOCKNUMBER = 15 & ELIGIBLE = 0 (as Palluck and AS did)

# all.schools.min <- all.schools[!(all.schools$OUTOFBLOCK.NUMERIC ==1 | all.schools$BLOCKNUMBER==15),]
all.schools.min <- all.schools[all.schools$ELIGIBLE==1,]


# Get numeric version of outcome (wearing orange band)
wrist.band.to.numeric <- str_extract(all.schools.min$WRISTOW2,"\\(?[0-9,.]+\\)?")
wrist.band.to.numeric <- as.numeric(gsub("\\((.+)\\)","\\1",wrist.band.to.numeric))
all.schools.min$WRISTOW2.NUMERIC <- wrist.band.to.numeric


# Create networks ---------------------------------------------------------

# PRE-PERIOD
# Create unique id
# schid*1000 + as.numeric(id)
all.schools.min$unique.id <- all.schools.min$SCHID*1000 + as.numeric(all.schools.min$ID)

attach(all.schools.min)
# Get questionnaire
survey.df <- as.data.frame(cbind(ST1, ST2, ST3, ST4, ST5, ST6, ST7, ST8, ST9, ST10))
survey.df <- apply(survey.df, 2, as.numeric)
# Get unique id
survey.df <- survey.df + SCHID*1000
# Save only eligible peers
survey.df.min <- apply(survey.df, 2,
                       function(x){ifelse(x %in% all.schools.min$unique.id,
                                          x,
                                          NA)})

# Get edges list
edgelist <- rbind(
  cbind(unique.id,survey.df.min[,1]),
  cbind(unique.id,survey.df.min[,2]),
  cbind(unique.id,survey.df.min[,3]),
  cbind(unique.id,survey.df.min[,4]),
  cbind(unique.id,survey.df.min[,5]),
  cbind(unique.id,survey.df.min[,6]),
  cbind(unique.id,survey.df.min[,7]),
  cbind(unique.id,survey.df.min[,8]),
  cbind(unique.id,survey.df.min[,9]),
  cbind(unique.id,survey.df.min[,10]),
  cbind(unique.id,unique.id)
)
# Remove NA
edgelist <- edgelist[!is.na(edgelist[,2]),]

either.peers.graph <- graph.data.frame(edgelist,directed=FALSE)
# final adj. matrix of `either` case (peers if at least one nominated the other as peer)
adj.matrix.either <- get.adjacency(either.peers.graph,sparse = FALSE)
adj.matrix.either[adj.matrix.either>1] <- 1
diag(adj.matrix.either) <- 0

### Best friends network
survey.bf.net <- as.data.frame(cbind(BF1, BF2))
survey.bf.net <- apply(survey.bf.net, 2, as.numeric)
# Get unique id
survey.bf.net <- survey.bf.net + SCHID*1000
# Save only eligible peers
survey.bf.net.min <- apply(survey.bf.net, 2,
                       function(x){ifelse(x %in% all.schools.min$unique.id,
                                          x,
                                          NA)})

# Get edges list
edgelist.bf.net <- rbind(
  cbind(unique.id,survey.bf.net.min[,1]),
  cbind(unique.id,survey.bf.net.min[,2]),
  cbind(unique.id,unique.id)
)

edgelist.bf.net <- edgelist.bf.net[!is.na(edgelist.bf.net[,2]),]

bf.either.peers.graph <- graph.data.frame(edgelist.bf.net,directed=FALSE)
bf.adj.matrix.either <- get.adjacency(bf.either.peers.graph,sparse = FALSE)
bf.adj.matrix.either[bf.adj.matrix.either>1] <- 1
diag(bf.adj.matrix.either) <- 0


### POST-PERIOD (WAVE 2)

# Get questionnaire
survey.df.W2 <- as.data.frame(cbind(ST1W2, ST2W2, ST3W2, ST4W2, ST5W2, ST6W2, ST7W2, ST8W2, ST9W2, ST10W2))
survey.df.W2 <- apply(survey.df.W2, 2, as.numeric)
# Get unique id
survey.df.W2 <- survey.df.W2 + SCHID*1000
# Save only eligible peers
survey.df.min.W2 <- apply(survey.df.W2, 2,
                       function(x){ifelse(x %in% all.schools.min$unique.id,
                                          x,
                                          NA)})

# Get edges list
edgelist.W2 <- rbind(
  cbind(unique.id,survey.df.min.W2[,1]),
  cbind(unique.id,survey.df.min.W2[,2]),
  cbind(unique.id,survey.df.min.W2[,3]),
  cbind(unique.id,survey.df.min.W2[,4]),
  cbind(unique.id,survey.df.min.W2[,5]),
  cbind(unique.id,survey.df.min.W2[,6]),
  cbind(unique.id,survey.df.min.W2[,7]),
  cbind(unique.id,survey.df.min.W2[,8]),
  cbind(unique.id,survey.df.min.W2[,9]),
  cbind(unique.id,survey.df.min.W2[,10]),
  cbind(unique.id,unique.id)
)

edgelist.W2 <- edgelist.W2[!is.na(edgelist.W2[,2]),]

either.peers.graph.W2 <- graph.data.frame(edgelist.W2,directed=FALSE)
# final adj. matrix of `either` case (peers if at least one nominated the other as peer)
adj.matrix.either.W2 <- get.adjacency(either.peers.graph.W2,sparse = FALSE)
adj.matrix.either.W2[adj.matrix.either.W2>1] <- 1
diag(adj.matrix.either.W2) <- 0


### Best friends network (W2)

survey.bf.net.W2 <- as.data.frame(cbind(BF1W2, BF2W2))
survey.bf.net.W2 <- apply(survey.bf.net.W2, 2, as.numeric)
# Get unique id
survey.bf.net.W2 <- survey.bf.net.W2 + SCHID*1000
# Save only eligible peers
survey.bf.net.W2.min <- apply(survey.bf.net.W2, 2,
                           function(x){ifelse(x %in% all.schools.min$unique.id,
                                              x,
                                              NA)})

# Get edges list
edgelist.bf.net.W2 <- rbind(
  cbind(unique.id,survey.bf.net.W2.min[,1]),
  cbind(unique.id,survey.bf.net.W2.min[,2]),
  cbind(unique.id,unique.id)
)

edgelist.bf.net.W2 <- edgelist.bf.net.W2[!is.na(edgelist.bf.net.W2[,2]),]

bf.either.peers.graph.W2 <- graph.data.frame(edgelist.bf.net.W2,directed=FALSE)
# final adj. matrix of `either` case (peers if at least one nominated the other as peer)
bf.adj.matrix.either.W2 <- get.adjacency(bf.either.peers.graph.W2,sparse = FALSE)
bf.adj.matrix.either.W2[bf.adj.matrix.either.W2>1] <- 1
diag(bf.adj.matrix.either.W2) <- 0


# Jaccard index -----------------------------------------------------------

net.list <- list(ST.network = adj.matrix.either, ST.W2.network = adj.matrix.either.W2,
                 BF.network = bf.adj.matrix.either, BF.W2.network = bf.adj.matrix.either.W2)

jacc.matrix <- matrix(NA,nrow=length(net.list),ncol=length(net.list))

for (i in seq(length(net.list))) {
  for (j in seq(length(net.list))) {
    jacc.matrix[i,j] <- jaccard_edgeset_similarity(A1 = net.list[[i]],
                                                   A2 = net.list[[j]])
  }
}

jacc.matrix[upper.tri(jacc.matrix)] <- NA
jacc.matrix <- round(jacc.matrix,3)
jacc.df <- as.data.frame(jacc.matrix)
colnames(jacc.df) <- names(net.list)
rownames(jacc.df) <- names(net.list)
options(knitr.kable.NA="")
kable(jacc.df, format = "latex", booktabs = T, digits = 3)

# Save networks data ------------------------------------------------------

saveRDS(list(ST.network = adj.matrix.either, ST.W2.network = adj.matrix.either.W2,
             BF.network = bf.adj.matrix.either, BF.W2.network = bf.adj.matrix.either.W2),
        "Reproducibility/Data_analyses/Palluck_et_al/adj_matrix_data.RDS")

# Estimate prob. matrices -------------------------------------------------

########
# This part can take a some time
########


expos.vec <- c("c111","c101","c011","c001","c000")
expos.contrast <- list(c("c111","c000"), c("c101","c000"),
                       c("c011","c000"), c("c001","c000"))
n <- nrow(adj.matrix.either)

set.seed(595957)

prob.mat.adj.either <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(adj.matrix.either),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)

prob.mat.adj.either.W2 <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(adj.matrix.either.W2),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)

prob.mat.adj.BF <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(bf.adj.matrix.either),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)

prob.mat.adj.BF.W2 <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(bf.adj.matrix.either.W2),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)

prob.mat.adj.either.and.BF <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(adj.matrix.either, bf.adj.matrix.either),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)

prob.mat.adj.either.and.BF.W2 <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(adj.matrix.either.W2, bf.adj.matrix.either.W2),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)


prob.mat.adj.either.W1.and.W2 <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(adj.matrix.either, adj.matrix.either.W2),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)

prob.mat.adj.all <- Get_prob_matrices_list(R = 10^4,
                                   n = n,
                                   Pz_function = treat.randomization.palluck,
                                   pz_func_args = list(n=n, pz=0.5, schid = SCHID),
                                   A.list = list(adj.matrix.either, adj.matrix.either.W2,
                                                 bf.adj.matrix.either, bf.adj.matrix.either.W2),
                                   exposures_contrast = expos.contrast,
                                   exposures_vec = expos.vec,
                                   threshold = NULL,
                                   Palluck.et.al = TRUE)
#
# saveRDS(list(prob.mat.adj.either = prob.mat.adj.either,
#              prob.mat.adj.either.W2 = prob.mat.adj.either.W2,
#              prob.mat.adj.BF = prob.mat.adj.BF,
#              prob.mat.adj.BF.W2 = prob.mat.adj.BF.W2,
#              prob.mat.adj.either.and.BF = prob.mat.adj.either.and.BF,
#              prob.mat.adj.either.W1.and.W2 = prob.mat.adj.either.W1.and.W2,
#              prob.mat.adj.all = prob.mat.adj.all),
#         "Simulations/data_analysis/Palluck_et_al/prob_mat_list.RDS")
# set.seed(985959)
# prob.mat.BF1.and.BF2 <- Get_prob_matrices_list(R = 10^4,
#                                    n = n,
#                                    Pz_function = treat.randomization.palluck,
#                                    pz_func_args = list(n=n, pz=0.5, schid = SCHID),
#                                    A.list = list(bf.adj.matrix.either, bf.adj.matrix.either.W2),
#                                    exposures_contrast = expos.contrast,
#                                    exposures_vec = expos.vec,
#                                    threshold = NULL,
#                                    Palluck.et.al = TRUE)
#
# saveRDS(list(prob.mat.BF1.and.BF2 = prob.mat.BF1.and.BF2),
#         "Simulations/data_analysis/Palluck_et_al/prob_mat_BF1.and.BF2.RDS")

# Estimate causal effects -------------------------------------------------

prob.mat <- readRDS("Reproducibility/Data_analyses/Palluck_et_al/prob_mat_list.RDS")
prob.mat.BFs <- readRDS("Reproducibility/Data_analyses/Palluck_et_al/prob_mat_BF1.and.BF2.RDS")

z.obs <- TREAT.NUMERIC
z.obs[z.obs==2] <- 0

# Compute CE for each of the adj. matrices list.

adj.either.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(adj.matrix.either),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat$prob.mat.adj.either,
                                       threshold = NULL,
                                       estimate.n = FALSE)

adj.either.W2.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(adj.matrix.either.W2),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat$prob.mat.adj.either.W2,
                                       threshold = NULL,
                                       estimate.n = FALSE)

adj.BF.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(bf.adj.matrix.either),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat$prob.mat.adj.BF,
                                       threshold = NULL,
                                       estimate.n = FALSE)

adj.BF.W2.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(bf.adj.matrix.either.W2),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat$prob.mat.adj.BF.W2,
                                       threshold = NULL,
                                       estimate.n = FALSE)

adj.either.and.BF.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(adj.matrix.either, bf.adj.matrix.either),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat$prob.mat.adj.either.and.BF,
                                       threshold = NULL,
                                       estimate.n = FALSE)

adj.either.W1.and.W2.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(adj.matrix.either, adj.matrix.either.W2),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat$prob.mat.adj.either.W1.and.W2,
                                       threshold = NULL,
                                       estimate.n = FALSE)

adj.BF1.and.BF2.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(bf.adj.matrix.either, bf.adj.matrix.either.W2),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat.BFs$prob.mat.BF1.and.BF2,
                                       threshold = NULL,
                                       estimate.n = FALSE)

adj.all.ce <- MR_CE_estimator_palluck(Z.obs = z.obs,
                                       school.treated = SCHTREAT.NUMERIC,
                                       Y.obs = WRISTOW2.NUMERIC,
                                       A.list = list(adj.matrix.either, adj.matrix.either.W2,
                                                      bf.adj.matrix.either, bf.adj.matrix.either.W2),
                                       exposures_contrast = expos.contrast,
                                       exposures_vec = expos.vec,
                                       Prob_matrices_list = prob.mat$prob.mat.adj.all,
                                       threshold = NULL,
                                       estimate.n = FALSE)

# combine all
c111.c000 <- rbindlist(list(adj.either.ce$`c111-c000`,
                              adj.either.W2.ce$`c111-c000`,
                              adj.BF.ce$`c111-c000`,
                              adj.BF.W2.ce$`c111-c000`,
                              adj.either.W1.and.W2.ce$`c111-c000`,
                              adj.either.and.BF.ce$`c111-c000`,
                              adj.BF1.and.BF2.ce$`c111-c000`,
                              adj.all.ce$`c111-c000`))
c111.c000$estimand <- c("c111-c000")
c111.c000$adj.mat <- c("ST","ST.W2","BF","BF.W2","ST.W1.and.W2",
                       "ST.and.BF", "BF.W1.and.BF.W2" ,"All")

c101.c000 <- rbindlist(list(adj.either.ce$`c101-c000`,
                              adj.either.W2.ce$`c101-c000`,
                              adj.BF.ce$`c101-c000`,
                              adj.BF.W2.ce$`c101-c000`,
                              adj.either.W1.and.W2.ce$`c101-c000`,
                              adj.either.and.BF.ce$`c101-c000`,
                              adj.BF1.and.BF2.ce$`c101-c000`,
                              adj.all.ce$`c101-c000`))
c101.c000$estimand <- c("c101-c000")
c101.c000$adj.mat <- c("ST","ST.W2","BF","BF.W2","ST.W1.and.W2",
                       "ST.and.BF", "BF.W1.and.BF.W2" ,"All")

c011.c000 <- rbindlist(list(adj.either.ce$`c011-c000`,
                              adj.either.W2.ce$`c011-c000`,
                              adj.BF.ce$`c011-c000`,
                              adj.BF.W2.ce$`c011-c000`,
                              adj.either.W1.and.W2.ce$`c011-c000`,
                              adj.either.and.BF.ce$`c011-c000`,
                              adj.BF1.and.BF2.ce$`c011-c000`,
                              adj.all.ce$`c011-c000`))
c011.c000$estimand <- c("c011-c000")
c011.c000$adj.mat <- c("ST","ST.W2","BF","BF.W2","ST.W1.and.W2",
                       "ST.and.BF", "BF.W1.and.BF.W2" ,"All")


c001.c000 <- rbindlist(list(adj.either.ce$`c001-c000`,
                              adj.either.W2.ce$`c001-c000`,
                              adj.BF.ce$`c001-c000`,
                              adj.BF.W2.ce$`c001-c000`,
                              adj.either.W1.and.W2.ce$`c001-c000`,
                              adj.either.and.BF.ce$`c001-c000`,
                              adj.BF1.and.BF2.ce$`c001-c000`,
                              adj.all.ce$`c001-c000`))
c001.c000$estimand <- c("c001-c000")
c001.c000$adj.mat <- c("ST","ST.W2","BF","BF.W2","ST.W1.and.W2",
                       "ST.and.BF", "BF.W1.and.BF.W2" ,"All")



palluck.et.al.results <- rbindlist(list(c001.c000,c011.c000,c101.c000,c111.c000))
palluck.et.al.results[,`:=`(ht_point_ci = paste0(round(ht_ce,3),
                                           " [",
                                           round(ht_ce - 1.96*sqrt(var_ht_ce),3),
                                           ", ",
                                           round(ht_ce + 1.96*sqrt(var_ht_ce),3),
                                           "]"),
                            hajek_point_ci = paste0(round(hajek_ce,3),
                                                 " [",
                                                 round(hajek_ce - 1.96*sqrt(var_hajek_ce),3),
                                                 ", ",
                                                 round(hajek_ce + 1.96*sqrt(var_hajek_ce),3),
                                                 "]"))]

palluck.et.al.results[,`:=`(ht_point_se = paste0(round(ht_ce,3),
                                           " (",
                                           round(sqrt(var_ht_ce),3),
                                           ")"),
                            hajek_point_se = paste0(round(hajek_ce,3),
                                                 " (",
                                                 round(sqrt(var_hajek_ce),3),
                                                 ")"))]

hajek.casted <- dcast.data.table(palluck.et.al.results[estimand%in%c("c101-c000","c001-c000"),],
                     " adj.mat ~ estimand", value.var = "hajek_point_se")

ht.casted <- dcast.data.table(palluck.et.al.results[estimand%in%c("c101-c000","c001-c000"),],
                     " adj.mat ~ estimand", value.var = "ht_point_se")

casted.results <- cbind(ht.casted,hajek.casted[,2:3])

library(kableExtra)

network_order <- c("ST","BF","ST.and.BF",
                   "ST.W2","BF.W2",
                   "ST.W1.and.W2",
                   "BF.W1.and.BF.W2",
                   "All")

network_labels <- c("ST","BF","ST & BF",
                   "ST 2","BF 2",
                   "ST 1 & 2",
                   "BF 1 & 2",
                   "ALL")

kable(x = casted.results[network_order,c(1,2,4,3,5)],
      format = "latex",
      booktabs = T,
      col.names = c("Networks",rep(c("HT", "Hajek"),2))) %>%
      # add_header_above(c(" " = 1, TeX("$\\tau(c001,c000)$") = 2,
      #                    TeX("$\\tau(c101,c000)$") = 2))
      add_header_above(c(" " = 1, "one" = 2,
                         "two" = 2))

network_order_plot <- c("ST","BF","ST.and.BF", "All")

network_labels_plot <- c("ST (pre)","BF (pre)","ST & BF (pre)", "ALL")

ht.results <- palluck.et.al.results[,.SD,.SDcols = c("ht_ce","var_ht_ce","estimand","adj.mat")]
ht.results[, method := "HT"]

hajek.results <- palluck.et.al.results[,.SD,.SDcols = c("hajek_ce","var_hajek_ce","estimand","adj.mat")]
hajek.results[, method := "Hajek"]

melted.results <- rbindlist(list(ht.results,hajek.results), use.names = FALSE)


library(ggplot2)
library(ggpubr)
library(latex2exp)


# ggplot(melted.results[estimand%in%c("c101-c000","c001-c000")],
ggplot(melted.results[method == "Hajek" & 
                        adj.mat %in% c("ST","BF","ST.and.BF", "All") &
                        estimand %in% c("c011-c000","c111-c000")],
       aes(x=factor(adj.mat,levels = rev(network_order_plot),
                    labels = rev(network_labels_plot)),
                           y = ht_ce,
           col = method,
           group = method,
           shape = method)) +
  geom_errorbar(aes(ymin = ht_ce - 1.96*sqrt(var_ht_ce),
                    ymax = ht_ce + 1.96*sqrt(var_ht_ce)),
                width = .4, position = position_dodge(0.3),
                linewidth =1.1, alpha = 0.7) +
  geom_point(position = position_dodge(0.3),
             cex = 8) +
  geom_hline(yintercept = 0, lty = "dashed",
             alpha = 0.6, linewidth = 1) +
  # scale_color_manual(values = c("HT" = "#990000","Hajek" = "#0065A9")) +
  scale_color_manual(values = c("Hajek" = "#0065A9")) +
  # scale_shape_manual(values = c("HT" = 16,"Hajek" = 18)) +
  labs(x="", y="") +
  scale_y_continuous(breaks = seq(-1,1,0.25),labels = seq(-1,1,0.25)) +
  coord_flip() +
  facet_wrap(~estimand, nrow = 1
             ,
             labeller = labeller(estimand = c("c111-c000" = "Overall Effect",
                                   "c011-c000" = "Indirect Effect"))
             ) +
  theme_pubclean() +
  theme(legend.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 22),
        axis.text.x = element_text(size = 18, face = "bold"),
        # legend.text = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 20),
        # legend.key.size = unit(1,"cm"),
        legend.position = "none",
        )


# Save figure
ggsave("Reproducibility/Data_analyses/Palluck_et_al/figures/Paluck_hajek_results.jpeg",
       width = 12, height = 6, units = "in")

# Appendix results --------------------------------------------------------


# hajek.casted.apdx <- dcast.data.table(palluck.et.al.results[!estimand%in%c("c101-c000","c001-c000"),],
#                                  " adj.mat ~ estimand", value.var = "hajek_point_se")
#
# ht.casted.apdx <- dcast.data.table(palluck.et.al.results[!estimand%in%c("c101-c000","c001-c000"),],
#                               " adj.mat ~ estimand", value.var = "ht_point_se")

hajek.casted.apdx <- dcast.data.table(palluck.et.al.results,
                                 # " adj.mat ~ estimand", value.var = "hajek_point_se")
                                 " adj.mat ~ estimand", value.var = "hajek_point_ci")

ht.casted.apdx <- dcast.data.table(palluck.et.al.results,
                              # " adj.mat ~ estimand", value.var = "ht_point_se")
                              " adj.mat ~ estimand", value.var = "ht_point_ci")

casted.results.apdx <- cbind(ht.casted.apdx,hajek.casted.apdx[,2:5])


kable(x = casted.results.apdx[network_order,c(1,2,6,3,7,4,8,5,9)],
      format = "latex",
      booktabs = T,
      col.names = c("Networks",rep(c("HT", "Hajek"),4))) %>%
  # add_header_above(c(" " = 1, TeX("$\\tau(c001,c000)$") = 2,
  #                    TeX("$\\tau(c101,c000)$") = 2))
  add_header_above(c(" " = 1, "one" = 2,
                     "two" = 2,
                     "three" = 2,
                     "four" = 2))


