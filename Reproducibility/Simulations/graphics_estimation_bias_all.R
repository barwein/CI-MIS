###
# Script that generate the graphic results of the bias simulation study
###

# Load libraries ----------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(stringi)
library(stringr)
library(kableExtra)


# Scenario I - Noise -----------------------------------------------------

# Read data

PA.noise.results <- fread("Reproducibility/Simulations/results/Bias_from_noise_PA_M1000.csv")

# Add true effects

PA.noise.results[,true_effect := rep(c(1,0.75,0.5,0.25,0.5),nrow(PA.noise.results)/5)]

# Melt and compute means

PA.noise.melted <- melt.data.table(PA.noise.results,
                                     id.vars = c("ce_contrast","param","true_effect"),
                                     measure.vars = c("ht_ce","hajek_ce","naive_diff_ce",
                                                      "exact.bias", "bounds.bias"))

PA.noise.summarized.wo.exact <- PA.noise.melted[!variable %in% c("exact.bias","bounds.bias"),
                                                .(mean_esti = abs(mean(value-true_effect)),
                                                  # perc025 = quantile(100*(value-true_effect)/true_effect,0.025),
                                                  # perc975 = quantile(100*(value-true_effect)/true_effect,0.975),
                                                  sd.h = sd(value-true_effect)),
                                               by = c("ce_contrast","param","variable")]

PA.noise.summarized.w.exact <- PA.noise.melted[variable %in% c("exact.bias","bounds.bias"),
                                                .(mean_esti = abs(mean(value)),
                                                  # perc025 = quantile(100*(value-true_effect)/true_effect,0.025),
                                                  # perc975 = quantile(100*(value-true_effect)/true_effect,0.975),
                                                  sd.h = sd(value)),
                                               by = c("ce_contrast","param","variable")]

PA.noise.summarized <- rbindlist(list(PA.noise.summarized.wo.exact, PA.noise.summarized.w.exact))


# PA.noise.summarized <- PA.noise.melted[,.(mean_esti = abs(mean(100*(value-true_effect)/true_effect)),
#                                       perc025 = quantile(100*(value-true_effect)/true_effect,0.025),
#                                       perc975 = quantile(100*(value-true_effect)/true_effect,0.975),
#                                       sd.h = sd(100*(value-true_effect)/true_effect)),
#                                    by = c("ce_contrast","param","variable")]
# 

PA.noise.plot.all <- ggplot(PA.noise.summarized, aes(x=factor(param),y=mean_esti, col = variable,
                                         group = variable, fill=variable)) +
  geom_point(position = position_dodge(0)) +
  geom_line(position = position_dodge(0)) +
  # geom_errorbar(aes(ymin = perc025, ymax = perc975),
  #               position = position_dodge(0.3),
  #               width = 0) +
  geom_hline(yintercept = 0, lty = "dashed") +
  facet_wrap("ce_contrast",nrow = 2) +
  theme_pubclean()

# Main text plot

PA.noise.plot.some <- ggplot(PA.noise.summarized[variable %in% c("ht_ce","hajek_ce") & 
                                                   # ce_contrast %in% c("c10-c00")],
                                                   ce_contrast %in% c("c11-c00")],
                       aes(x = factor(param),
                           y = mean_esti,
                           color = variable,
                           group = variable,
                           fill = variable,
                           shape = variable)) +
  geom_line(linewidth = 1, alpha = 0.7, show.legend = F) +
  geom_point(size = 12) +
  scale_fill_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9"),
                    labels = c("HT","Hajek")) +
  scale_color_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9"),
                     labels = c("HT","Hajek")) +
  scale_shape_manual(values = c("ht_ce" = 16,"hajek_ce" = 18),
                     labels = c("HT","Hajek")) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  # scale_y_continuous(breaks = seq(0,.1,0.025), limits = c(0,.1)) +
  scale_y_continuous(breaks = seq(0,.3,0.05), limits = c(0,.3)) +
  # scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20)) +
  # labs(x=TeX("$\\eta$"), y = "Bias (%)", title = "(I)") +
  labs(x=TeX("$\\eta$"), y = "Abs. bias", title = "(I)") +
  # facet_wrap(~ce_contrast,nrow = 1,
  #            labeller = labeller(ce_contrast = c("c10-c00"="Direct effect"))) +
  guides(fill = guide_legend(override.aes = list(size = 12))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size =24, face = "bold", angle = 30),
        axis.text.y = element_text(size =24, face = "bold"),
        axis.title.x = element_text(size = 44, face = "bold"),
        axis.title.y = element_text(size = 26, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=28, face = "bold"),
        legend.position = "bottom",
        # legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=28, face="bold"),
        legend.key.size = unit(1, "cm"),  # Increase key (symbol) size
        legend.box.spacing = unit(0.3, "cm"),  # Adjust the box spacing
        # legend.background = element_blank(),  # Remove the background
        legend.background = element_rect(fill = "white"),  # Remove the background
        # legend.position = c(0.5,0.1),
        # legend.text = element_text(size = 24, face = "bold"),
        # legend.key.size = unit(1.2,"cm"),
        plot.title = element_text(size = 36, face = "bold", hjust = .5, vjust = -0.5),
        plot.margin=unit(c(1,0.15,1,1),"cm"))

# Appendix plots

PA.noise.summarized$labels = factor(PA.noise.summarized$ce_contrast,
                                    labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))

PA.noise.plot.apdx <- ggplot(PA.noise.summarized[!variable %in% c("naive_diff_ce","bounds.bias")
                                                 ],
                                                 # & !ce_contrast %in% c("c10-c00")],
                       aes(x = factor(param),
                           y = mean_esti,
                           color = variable,
                           group = variable,
                           fill = variable,
                           shape = variable)) +
                  geom_line(linewidth = 1, alpha = 0.7, show.legend = F) +
                  geom_point(size = 7,alpha=0.6) + 
                  scale_fill_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9","exact.bias"="darkgreen"),
                                    labels = c("HT","Hajek","Exact bias")) +
                  scale_color_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9","exact.bias"="darkgreen"),
                                     labels = c("HT","Hajek","Exact bias")) +
                  scale_shape_manual(values = c("ht_ce" = 16,"hajek_ce" = 18,"exact.bias"=17),
                                     labels = c("HT","Hajek","Exact bias")) +
                  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
                  # scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20)) +
                  # labs(x=TeX("$\\eta$"), y = "Bias (%)") +
                  labs(x=TeX("$\\eta$"), y = "Abs. bias") +
                  facet_wrap(~labels,nrow = 1, labeller = label_parsed, scales = "free_y") +
                  guides(fill = guide_legend(override.aes = list(size = 10))) +
                  theme_pubclean() +
                  theme(axis.text.x = element_text(size =20, face = "bold",angle=-90,hjust = 0.8),
                        axis.text.y = element_text(size =20, face = "bold"),
                        axis.title.x = element_text(size = 38, face = "bold"),
                        axis.title.y = element_text(size = 20, face="bold"),
                        strip.background = element_blank(),
                        strip.text = element_text(size=24, face = "bold"),
                        legend.position = "top",
                        # legend.position = "none",
                        # legend.position = c(0.5,0.1),
                        legend.title = element_blank(),
                        legend.text = element_text(size = 24, face = "bold"),
                        legend.key.size = unit(1.2,"cm"),
                        # plot.title = element_text(size = 34, face = "bold", hjust = .5, vjust = -0.5)
                        )

ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_PA_noise_bias_1000iter.jpeg",
       plot = PA.noise.plot.apdx, 
       width = 18, height = 10)


# Number of misclassified exposures

n.c11 <- str_split(str_extract_all(string = PA.noise.results$n.expos.misclass, pattern = "m.c11=\\d+"),
                   "=")
n.c11 <- as.numeric(sapply(n.c11,"[",2))

n.c01 <- str_split(str_extract_all(string = PA.noise.results$n.expos.misclass, pattern = "m.c01=\\d+"),
                   "=")
n.c01 <- as.numeric(sapply(n.c01,"[",2))

n.c10 <- str_split(str_extract_all(string = PA.noise.results$n.expos.misclass, pattern = "m.c10=\\d+"),
                   "=")
n.c10 <- as.numeric(sapply(n.c10,"[",2))

n.c00 <- str_split(str_extract_all(string = PA.noise.results$n.expos.misclass, pattern = "m.c00=\\d+"),
                   "=")
n.c00 <- as.numeric(sapply(n.c00,"[",2))

PA.noise.results[,`:=`(n.c11 = n.c11, n.c01 = n.c01, n.c10 = n.c10, n.c00 = n.c00)]

PA_noise_n_miclass_melted <- melt.data.table(PA.noise.results[ce_contrast=="c11-c00",],
                                       id.vars = "param",
                                       measure.vars = c("n.c11","n.c01","n.c10","n.c00"))

PA_noise_n_miclass_mean <- PA.noise.results[ce_contrast=="c11-c00",
                                .(mean.c11 = mean(n.c11),
                                  mean.c01 = mean(n.c01),
                                  mean.c10 = mean(n.c10),
                                  mean.c00 = mean(n.c00)),
                                by = "param"]

expos.labels <- c("c11","c01","c10","c00")
names(expos.labels) <- c("n.c11","n.c01","n.c10","n.c00")


PA_noise_n_miclass_melted$labels = factor(PA_noise_n_miclass_melted$variable,
                                          labels = c(TeX(r"($c_{\1\1}$)",output = "character"),
                                                     TeX(r"($c_{\0\1}$)",output = "character"),
                                                     TeX(r"($c_{\1\0}$)",output = "character"),
                                                     TeX(r"($c_{\0\0}$)",output = "character")))

PA.noise.n.misclass.plot <- ggplot(PA_noise_n_miclass_melted,
                             aes(x = factor(param),
                                 y = value,
                                 # color = variable,
                                 group = param,
                                 # fill = variable,
                                 # shape = variable
                             )) +
  geom_boxplot(fill="darkgrey", alpha = 0.7) +
  labs(x=TeX("$\\eta$"), y = "# Missclassified exposures") +
  facet_wrap(~labels, nrow = 2,
             labeller = label_parsed) +
  # theme_pubclean() +
  theme_bw() +
  theme(axis.text.x = element_text(size =18, face = "bold"),
        axis.text.y = element_text(size =18, face = "bold"),
        axis.title.x = element_text(size = 35, face = "bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=24, face = "bold"))

ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_PA_noise_n_missclassified_exposures.jpeg",
       plot = PA.noise.n.misclass.plot, 
       width = 18, height = 10)

# HT and Hajek empirical SD

PA_noise_results.sd <- PA.noise.results[,.(mean.ht = mean(ht_ce-true_effect),
                               sd.ht = sd(ht_ce-true_effect),
                               mean.hajek = mean(hajek_ce-true_effect),
                               sd.hajek = sd(hajek_ce-true_effect)),
                            by = c("ce_contrast","param")]


PA_noise_results.sd <- melt.data.table(PA_noise_results.sd,
                                id.vars = c("ce_contrast","param"),
                                measure.vars = c("sd.ht","sd.hajek"))

PA_noise_results.sd$labels = factor(PA_noise_results.sd$ce_contrast,
                                    labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))

PA.noise.sd.plot <- ggplot(PA_noise_results.sd,
                     aes(x = factor(param),
                         y = value,
                         color = variable,
                         group = variable,
                         fill = variable,
                         shape = variable)) +
  geom_line(linewidth = 1, alpha = 0.7, show.legend = F) +
  geom_point(size = 7, alpha=0.6) +
  scale_fill_manual(values = c("sd.ht" = "#990000","sd.hajek" = "#0065A9"),
                    labels = c("HT","Hajek")) +
  scale_color_manual(values = c("sd.ht" = "#990000","sd.hajek" = "#0065A9"),
                     labels = c("HT","Hajek")) +
  scale_shape_manual(values = c("sd.ht" = 16,"sd.hajek" = 18),
                     labels = c("HT","Hajek")) +
  # geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  labs(x=TeX("$\\eta$"), y = "Empirical SD") +
  # facet_wrap(~ce_contrast,nrow = 1,
  #            labeller = labeller(ce_contrast = c("c01-c00"="Indirect effect",
  #                                                "c11-c00"="Total effect"))) +
  facet_wrap(~labels, nrow = 1, labeller = label_parsed) +
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  # theme_pubclean() +
  theme_bw() +
  theme(axis.text.x = element_text(size =18, face = "bold", angle = -90),
        axis.text.y = element_text(size =18, face = "bold"),
        axis.title.x = element_text(size = 35, face = "bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=24, face = "bold"),
        legend.position = "top",
        # legend.position = c(0.5,0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1.2,"cm"))

ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_PA_noise_SD_plot.jpeg",
       plot = PA.noise.sd.plot, 
       width = 18, height = 10)

# Network similarity with Jaccard index

PA.noise.jaccard <- rbind(unique(PA.noise.results$param),
                          round(unique(PA.noise.results$jac),3))
rownames(PA.noise.jaccard) <- c("eta","J")

kable(x = PA.noise.jaccard,format = "latex", booktabs = T )


# Spaghetti plots


# csv_names <- dir("Reproducibility/Simulations/results/APDX_noise/")
# path_name <- "Reproducibility/Simulations/results/APDX_noise/"
# 
# combined_dt <- data.table()
# 
# for (i in seq(length(csv_names))) {
#   
#   res_name <- csv_names[i]
#   curr_dt <- fread(paste0(path_name,res_name))
#   curr_dt$m.it <- i
#   combined_dt <- rbindlist(list(combined_dt,curr_dt))
# }
# 
# 
# write.csv(combined_dt,
#           "Reproducibility/Simulations/results/APDX_noise_spagetti.csv",
#           row.names = FALSE)


spaghetti.noise <- fread("Reproducibility/Simulations/results/APDX_noise_spagetti.csv")

spaghetti.noise[,true_effect := rep(c(1,0.75,0.5,0.25,0.5),nrow(spaghetti.noise)/5)]

PA.noise.results[,m.it := "A"]

spaghetti.noise <- rbindlist(list(spaghetti.noise,
                                  PA.noise.results[,.SD,.SDcols = names(spaghetti.noise)]))


spaghetti.noise.melted <- melt.data.table(spaghetti.noise,
                                   id.vars = c("ce_contrast","param","true_effect","m.it"),
                                   measure.vars = c("ht_ce","hajek_ce"))


spaghetti.noise.summarized <- spaghetti.noise.melted[,
                                                .(mean_esti = abs(mean(value-true_effect)),
                                                  # perc025 = quantile(100*(value-true_effect)/true_effect,0.025),
                                                  # perc975 = quantile(100*(value-true_effect)/true_effect,0.975),
                                                  sd.h = sd(value-true_effect)),
                                                by = c("ce_contrast","param","variable","m.it")]

spaghetti.noise.summarized[,apdx := ifelse(m.it == "A",FALSE,TRUE)]

spaghetti.noise.summarized$labels = factor(spaghetti.noise.summarized$ce_contrast,
                                      labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                                 TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                                 TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                                 TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                                 TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))

noise.spagehtti.plot <- ggplot(spaghetti.noise.summarized[ce_contrast %in% c("c11-c00","c11-c10") &
                                                       variable == "hajek_ce",],
                               aes(x = factor(param), 
                                   y = mean_esti,
                                   color = apdx,
                                   group = m.it,
                                   fill = apdx, 
                                   alpha = apdx,
                                   linewidth = apdx)) +
  geom_line(show.legend = F) + 
  scale_alpha_manual(values = c(0.8,0.2)) +
  scale_linewidth_discrete(range = c(1.6,1)) +
  # geom_point(size = 12, shape = 18) + 
  scale_fill_manual(values = c("TRUE" = "gray44","FALSE" = "#0065A9")) +
  # ,labels = c("HT","Hajek")) +
  scale_color_manual(values = c("TRUE" = "gray44","FALSE" = "#0065A9")) +
  # ,labels = c("HT","Hajek")) +
  # scale_shape_manual(values = c("ht_ce" = 16,"hajek_ce" = 18),
  #                    labels = c("HT","Hajek")) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  scale_y_continuous(breaks = seq(0,0.4,0.1),limits = c(0,0.4)) +
  labs(x=TeX("$\\eta$"), y = "Abs. bias") +
  facet_wrap(~labels,nrow = 1, labeller = label_parsed) +
  guides(fill = guide_legend(override.aes = list(size = 12))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size =26, face = "bold"),
        axis.text.y = element_text(size =26, face = "bold"),
        axis.title.x = element_text(size = 44, face = "bold"),
        axis.title.y = element_text(size = 26, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=28, face = "bold"),
        # legend.position = "top",
        legend.position = "none",
        # legend.position = c(0.5,0.1),
        # legend.title = element_blank(),
        # legend.text = element_text(size = 24, face = "bold"),
        # legend.key.size = unit(1.2,"cm"),
        plot.title = element_text(size = 36, face = "bold", hjust = .5, vjust = -0.5))


ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_noise_spagehti.jpeg",
       plot = noise.spagehtti.plot, 
       width = 20, height = 12)


# Scenario II - Censoring ---------------------------------------------------------------

# Read data

PA.censor.results <- fread("Reproducibility/Simulations/results/Bias_from_censoring_PA_M1000.csv")

# Add true effects

PA.censor.results[,true_effect := rep(c(1,0.75,0.5,0.25,0.5),nrow(PA.censor.results)/5)]

# Melt and compute means

PA.censor.melted <- melt.data.table(PA.censor.results,
                                   id.vars = c("ce_contrast","param","true_effect"),
                                   measure.vars = c("ht_ce","hajek_ce","naive_diff_ce",
                                                    "exact.bias","bounds.bias"))

PA.censor.summarized.wo.exact <- PA.censor.melted[!variable %in% c("exact.bias","bounds.bias")
                                          ,.(mean_esti = abs(mean(value-true_effect)),
                                             sd.h = sd(value-true_effect)),
                                       by = c("ce_contrast","param","variable")]

PA.censor.summarized.w.exact <- PA.censor.melted[variable %in% c("exact.bias","bounds.bias")
                                          ,.(mean_esti = abs(mean(value)),
                                             sd.h = sd(value)),
                                       by = c("ce_contrast","param","variable")]

PA.censor.summarized <- rbindlist(list(PA.censor.summarized.wo.exact, PA.censor.summarized.w.exact))


# PA.censor.summarized <- PA.censor.melted[,.(mean_esti = abs(mean(100*(value-true_effect)/true_effect)),
#                                           perc025 = quantile(100*(value-true_effect)/true_effect,0.025),
#                                           perc975 = quantile(100*(value-true_effect)/true_effect,0.975),
#                                           sd.h = sd(100*(value-true_effect)/true_effect)),
#                                        by = c("ce_contrast","param","variable")]
# 

PA.censor.plot.all <- ggplot(PA.censor.summarized, aes(x=factor(param,levels = as.character(seq(7,1))),
                                                       y=mean_esti, col = variable,
                                                     group = variable, fill=variable)) +
  geom_point(position = position_dodge(0)) +
  geom_line(position = position_dodge(0)) +
  # geom_errorbar(aes(ymin = perc025, ymax = perc975),
  #               position = position_dodge(0.3),
  #               width = 0) +
  geom_hline(yintercept = 0, lty = "dashed") +
  facet_wrap("ce_contrast",nrow = 2,scales = "free_y") +
  theme_pubclean()

# Main text plot

PA.censor.plot.some <- ggplot(PA.censor.summarized[variable %in% c("ht_ce","hajek_ce") & 
                                                     ce_contrast %in% c("c10-c00")],
                             aes(x = factor(param,levels = as.character(seq(7,1))),
                                 y = mean_esti,
                                 color = variable,
                                 group = variable,
                                 fill = variable,
                                 shape = variable)) +
  geom_line(linewidth = 1, alpha = 0.7, show.legend = F) +
  geom_point(size = 12) +
  scale_fill_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9"),
                    labels = c("HT","Hajek")) +
  scale_color_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9"),
                     labels = c("HT","Hajek")) +
  scale_shape_manual(values = c("ht_ce" = 16,"hajek_ce" = 18),
                     labels = c("HT","Hajek")) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  scale_y_continuous(breaks = seq(0,.1,0.025), limits = c(0,.1)) +
  # scale_y_continuous(breaks = seq(0,20,5), limits = c(0,20)) +
  # labs(x="K", y = "Bias (%)", title = "(II)") +
  labs(x="K", y = "Abs. bias", title = "(II)") +
  # facet_wrap(~ce_contrast,nrow = 1,
  #            labeller = labeller(ce_contrast = c("c10-c00"="Direct effect"))) +
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size =24, face = "bold"),
        axis.text.y = element_text(size =24, face = "bold"),
        # axis.text.y = element_blank(),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 26, face="bold"),
        # axis.title.y = element_text(size = 22, face="bold"),
        # axis.title.y = element_blank(),
        # axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=28, face = "bold"),
        legend.position = "bottom",
        # legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=28, face="bold"),
        legend.key.size = unit(1, "cm"),  # Increase key (symbol) size
        legend.box.spacing = unit(0.3, "cm"),  # Adjust the box spacing
        # legend.background = element_blank(),  # Remove the background
        legend.background = element_rect(fill = "white"),  # Remove the background
        # legend.position = c(0.5,0.1),
        # legend.title = element_blank(),
        # legend.text = element_text(size = 24, face = "bold"),
        # legend.key.size = unit(1.2,"cm"),
        plot.title = element_text(size = 36, face = "bold", hjust = .5, vjust = -0.5),
        plot.margin=unit(c(1,1,1,0.15),"cm"))



plot.noise.censored.combined <- ggarrange(PA.noise.plot.some,PA.censor.plot.some,
                                ncol = 2, nrow = 1, align = "h",
                                # widths = c(1.57,1)
                                widths = c(1.25,1),
                                # widths = c(1.75,1)
                                common.legend = TRUE, legend = "bottom"
                                )

ggsave(filename = "Reproducibility/Simulations/graphics/main/Bias_noise_censor_combined.jpeg",
       plot = plot.noise.censored.combined, 
       dpi = 500,
       width = 18, height = 10)


# Appendix plots

PA.censor.summarized$labels = factor(PA.censor.summarized$ce_contrast,
                                     labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))

PA.censor.plot.apdx <- ggplot(PA.censor.summarized[variable %in% c("ht_ce","hajek_ce","exact.bias")],
                             aes(x = factor(param,levels = as.character(seq(7,1))),
                                 y = mean_esti,
                                 color = variable,
                                 group = variable,
                                 fill = variable,
                                 shape = variable)) +
  geom_line(linewidth = 1, alpha = 0.7, show.legend = F) +
  geom_point(size = 7, alpha=0.6) +
  scale_fill_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9","exact.bias"="darkgreen"),
                    labels = c("HT","Hajek","Exact bias")) +
  scale_color_manual(values = c("ht_ce" = "#990000","hajek_ce" = "#0065A9","exact.bias"="darkgreen"),
                     labels = c("HT","Hajek","Exact bias")) +
  scale_shape_manual(values = c("ht_ce" = 16,"hajek_ce" = 18,"exact.bias"=17),
                     labels = c("HT","Hajek","Exact bias")) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  # scale_y_continuous(breaks = seq(0,280,20), limits = c(0,300)) +
  labs(x=TeX("$\\eta$"), y = "Abs. bias") +
  facet_wrap(~labels,nrow = 1,labeller = label_parsed, scales = "free_y") +
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size =20, face = "bold"),
        axis.text.y = element_text(size =20, face = "bold"),
        axis.title.x = element_text(size = 38, face = "bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=24, face = "bold"),
        legend.position = "top",
        # legend.position = "none",
        # legend.position = c(0.5,0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, face = "bold"),
        legend.key.size = unit(1.2,"cm"),
        # plot.title = element_text(size = 34, face = "bold", hjust = .5, vjust = -0.5)
  )

ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_PA_censor_bias_1000iter.jpeg",
       plot = PA.censor.plot.apdx, 
       width = 18, height = 10)


# Number of misclassified exposures

n.c11 <- str_split(str_extract_all(string = PA.censor.results$n.expos.misclass, pattern = "m.c11=\\d+"),
                   "=")
n.c11 <- as.numeric(sapply(n.c11,"[",2))

n.c01 <- str_split(str_extract_all(string = PA.censor.results$n.expos.misclass, pattern = "m.c01=\\d+"),
                   "=")
n.c01 <- as.numeric(sapply(n.c01,"[",2))

n.c10 <- str_split(str_extract_all(string = PA.censor.results$n.expos.misclass, pattern = "m.c10=\\d+"),
                   "=")
n.c10 <- as.numeric(sapply(n.c10,"[",2))

n.c00 <- str_split(str_extract_all(string = PA.censor.results$n.expos.misclass, pattern = "m.c00=\\d+"),
                   "=")
n.c00 <- as.numeric(sapply(n.c00,"[",2))

PA.censor.results[,`:=`(n.c11 = n.c11, n.c01 = n.c01, n.c10 = n.c10, n.c00 = n.c00)]

PA_censor_n_miclass_melted <- melt.data.table(PA.censor.results[ce_contrast=="c11-c00",],
                                             id.vars = "param",
                                             measure.vars = c("n.c11","n.c01","n.c10","n.c00"))

PA_censor_n_miclass_mean <- PA.censor.results[ce_contrast=="c11-c00",
                                            .(mean.c11 = mean(n.c11),
                                              mean.c01 = mean(n.c01),
                                              mean.c10 = mean(n.c10),
                                              mean.c00 = mean(n.c00)),
                                            by = "param"]

expos.labels <- c("c11","c01","c10","c00")
names(expos.labels) <- c("n.c11","n.c01","n.c10","n.c00")

PA_censor_n_miclass_melted$labels = factor(PA_censor_n_miclass_melted$variable,
                                           labels = c(TeX(r"($c_{\1\1}$)",output = "character"),
                                                      TeX(r"($c_{\0\1}$)",output = "character"),
                                                      TeX(r"($c_{\1\0}$)",output = "character"),
                                                      TeX(r"($c_{\0\0}$)",output = "character")))

PA.censor.n.misclass.plot <- ggplot(PA_censor_n_miclass_melted,
                                   aes(x = factor(param,levels = as.character(seq(7,1))),
                                       y = value,
                                       # color = variable,
                                       group = param,
                                       # fill = variable,
                                       # shape = variable
                                   )) +
  geom_boxplot(fill="darkgrey", alpha = 0.7) +
  labs(x=TeX("$\\eta$"), y = "# Missclassified exposures") +
  facet_wrap(~labels, nrow = 2,
             labeller = label_parsed) +
  # theme_pubclean() +
  theme_bw() +
  theme(axis.text.x = element_text(size =18, face = "bold"),
        axis.text.y = element_text(size =18, face = "bold"),
        axis.title.x = element_text(size = 35, face = "bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=24, face = "bold"))

ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_PA_censor_n_missclassified_exposures.jpeg",
       plot = PA.censor.n.misclass.plot, 
       width = 18, height = 10)

# HT and Hajek empirical SD

PA_censor_results.sd <- PA.censor.results[,.(mean.ht = mean(ht_ce-true_effect),
                                           sd.ht = sd(ht_ce-true_effect),
                                           mean.hajek = mean(hajek_ce-true_effect),
                                           sd.hajek = sd(hajek_ce-true_effect)),
                                        by = c("ce_contrast","param")]


PA_censor_results.sd <- melt.data.table(PA_censor_results.sd,
                                       id.vars = c("ce_contrast","param"),
                                       measure.vars = c("sd.ht","sd.hajek"))


PA_censor_results.sd$labels = factor(PA_censor_results.sd$ce_contrast,
                                     labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                                TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))

PA.censor.sd.plot <- ggplot(PA_censor_results.sd,
                           aes(x = factor(param,levels = as.character(seq(7,1))),
                               y = value,
                               color = variable,
                               group = variable,
                               fill = variable,
                               shape = variable)) +
  geom_line(linewidth = 1, alpha = 0.7, show.legend = F) +
  geom_point(size = 7,alpha=0.6) +
  scale_fill_manual(values = c("sd.ht" = "#990000","sd.hajek" = "#0065A9"),
                    labels = c("HT","Hajek")) +
  scale_color_manual(values = c("sd.ht" = "#990000","sd.hajek" = "#0065A9"),
                     labels = c("HT","Hajek")) +
  scale_shape_manual(values = c("sd.ht" = 16,"sd.hajek" = 18),
                     labels = c("HT","Hajek")) +
  # geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  labs(x=TeX("$\\eta$"), y = "Empirical SD") +
  # facet_wrap(~ce_contrast,nrow = 1,
  #            labeller = labeller(ce_contrast = c("c01-c00"="Indirect effect",
  #                                                "c11-c00"="Total effect"))) +
  facet_wrap(~labels, nrow = 1, labeller = label_parsed) +
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  # theme_pubclean() +
  theme_bw() +
  theme(axis.text.x = element_text(size =18, face = "bold"),
        axis.text.y = element_text(size =18, face = "bold"),
        axis.title.x = element_text(size = 35, face = "bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=24, face = "bold"),
        legend.position = "top",
        # legend.position = c(0.5,0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, face = "bold"),
        legend.key.size = unit(1.2,"cm"))

ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_PA_censor_SD_plot.jpeg",
       plot = PA.censor.sd.plot, 
       width = 18, height = 10)

# Network similarity with Jaccard index

PA.censor.jaccard <- rbind(order(unique(PA.censor.results$param),decreasing = T),
                          round(unique(PA.censor.results$jac),3)[seq(7,1)])
rownames(PA.censor.jaccard) <- c("eta","J")

kable(x = PA.censor.jaccard,format = "latex", booktabs = T )


# Spaghetti plots


# csv_names <- dir("Reproducibility/Simulations/results/APDX_censor/")
# path_name <- "Reproducibility/Simulations/results/APDX_censor/"
# 
# combined_dt <- data.table()
# 
# for (i in seq(length(csv_names))) {
# 
#   res_name <- csv_names[i]
#   curr_dt <- fread(paste0(path_name,res_name))
#   curr_dt$m.it <- i
#   combined_dt <- rbindlist(list(combined_dt,curr_dt))
# }
# 
# 
# write.csv(combined_dt,
#           "Reproducibility/Simulations/results/APDX_censor_spagetti.csv",
#           row.names = FALSE)


spaghetti.censor <- fread("Reproducibility/Simulations/results/APDX_censor_spagetti.csv")

spaghetti.censor[,true_effect := rep(c(1,0.75,0.5,0.25,0.5),nrow(spaghetti.censor)/5)]

PA.censor.results[,m.it := "A"]

spaghetti.censor <- rbindlist(list(spaghetti.censor,
                                  PA.censor.results[,.SD,.SDcols = names(spaghetti.censor)]))


spaghetti.censor.melted <- melt.data.table(spaghetti.censor,
                                          id.vars = c("ce_contrast","param","true_effect","m.it"),
                                          measure.vars = c("ht_ce","hajek_ce"))


spaghetti.censor.summarized <- spaghetti.censor.melted[,
                                                     .(mean_esti = abs(mean(value-true_effect)),
                                                       # perc025 = quantile(100*(value-true_effect)/true_effect,0.025),
                                                       # perc975 = quantile(100*(value-true_effect)/true_effect,0.975),
                                                       sd.h = sd(value-true_effect)),
                                                     by = c("ce_contrast","param","variable","m.it")]

spaghetti.censor.summarized[,apdx := ifelse(m.it == "A",FALSE,TRUE)]

spaghetti.censor.summarized$labels = factor(spaghetti.censor.summarized$ce_contrast,
                                           labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                                      TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                                      TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                                      TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                                      TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))

censor.spagehtti.plot <- ggplot(spaghetti.censor.summarized[ce_contrast %in% c("c10-c00","c11-c00") &
                                                            variable == "hajek_ce",],
                               aes(x = factor(param,levels = as.character(seq(7,1))), 
                                   y = mean_esti,
                                   color = apdx,
                                   group = m.it,
                                   fill = apdx, 
                                   alpha = apdx,
                                   linewidth = apdx)) +
  geom_line(show.legend = F) + 
  scale_alpha_manual(values = c(0.8,0.2)) +
  scale_linewidth_discrete(range = c(1.6,1)) +
  # geom_point(size = 12, shape = 18) + 
  scale_fill_manual(values = c("TRUE" = "gray44","FALSE" = "#0065A9")) +
  # ,labels = c("HT","Hajek")) +
  scale_color_manual(values = c("TRUE" = "gray44","FALSE" = "#0065A9")) +
  # ,labels = c("HT","Hajek")) +
  # scale_shape_manual(values = c("ht_ce" = 16,"hajek_ce" = 18),
  #                    labels = c("HT","Hajek")) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  # scale_y_continuous(breaks = seq(0,0.4,0.1),limits = c(0,0.4)) +
  labs(x="K", y = "Abs. bias") +
  facet_wrap(~labels,nrow = 1, labeller = label_parsed) +
  guides(fill = guide_legend(override.aes = list(size = 12))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size =26, face = "bold"),
        axis.text.y = element_text(size =26, face = "bold"),
        axis.title.x = element_text(size = 44, face = "bold"),
        axis.title.y = element_text(size = 26, face="bold"),
        strip.background = element_blank(),
        strip.text = element_text(size=28, face = "bold"),
        # legend.position = "top",
        legend.position = "none",
        # legend.position = c(0.5,0.1),
        # legend.title = element_blank(),
        # legend.text = element_text(size = 24, face = "bold"),
        # legend.key.size = unit(1.2,"cm"),
        plot.title = element_text(size = 36, face = "bold", hjust = .5, vjust = -0.5))


ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/APDX_censor_spagehti.jpeg",
       plot = censor.spagehtti.plot, 
       width = 20, height = 12)




# Scenario III - Contamination -----------------------------------------------------------

# Read data

CR_fixed <- fread("Reproducibility/Simulations/results/CR_bias_Nc500_fixed_M1000.csv")
CR_vari <- fread("Reproducibility/Simulations/results/CR_bias_Nc2000_vari_M1000.csv")


# Summarized

cr_fixed_summarized <- CR_fixed[,.(mean_bias = abs(mean(y_hat-y_true))),
                                   # perc025 = quantile(100*(y_hat-y_true)/y_true,0.025),
                                   # perc975 = quantile(100*(y_hat-y_true)/y_true,0.975)),
                                by = "theta_"]

cr_vari_summarized <- CR_vari[,.(mean_bias = abs(mean(y_hat-y_true))),
                                 # perc025 = quantile(100*(y_hat-y_true)/y_true,0.025),
                                 # perc975 = quantile(100*(y_hat-y_true)/y_true,0.975)),
                              by = "theta_"]

# cr_fixed_summarized <- CR_fixed[,.(mean_bias = abs(mean(100*(y_hat-y_true)/y_true)),
#                                    perc025 = quantile(100*(y_hat-y_true)/y_true,0.025),
#                                    perc975 = quantile(100*(y_hat-y_true)/y_true,0.975)),
#                                 by = "theta_"]
# 
# cr_vari_summarized <- CR_vari[,.(mean_bias = abs(mean(100*(y_hat-y_true)/y_true)),
#                                  perc025 = quantile(100*(y_hat-y_true)/y_true,0.025),
#                                  perc975 = quantile(100*(y_hat-y_true)/y_true,0.975)),
#                               by = "theta_"]
# 
cr_fixed_summarized[,cluster_type := "Fixed"]
cr_vari_summarized[,cluster_type := "Varied"]

cr_summarized <- rbindlist(list(cr_fixed_summarized, cr_vari_summarized))


bias_cluster_plot <- 
  ggplot(cr_summarized[theta_ %in% seq(0,0.05,0.005),],
  # ggplot(cr_summarized[theta_ %in% seq(0,0.1,0.01),],
                            aes(x=theta_, y = mean_bias,
                                fill = cluster_type, color = cluster_type,
                                shape = cluster_type)) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 1.1) +
  # geom_point(size = 11) +
  geom_point(size = 10) +
  geom_line(linewidth = 1, alpha = 0.7, show.legend = F) +
  scale_x_continuous(breaks = seq(0,0.05,0.005), labels = as.character(seq(0,0.05,0.005))) +
  # scale_x_continuous(breaks = seq(0,0.1,0.01), labels = as.character(seq(0,0.1,0.01))) +
  # scale_y_continuous(breaks = seq(0,40,10), limits = c(0,40)) +
  scale_y_continuous(breaks = seq(0,0.35,0.05), limits = c(0,.35)) +
  scale_color_manual(values = c("Fixed" = "#009933", "Varied" = "#660033")) +
  scale_fill_manual(values = c("Fixed" = "#009933", "Varied" = "#660033")) +
  scale_shape_manual(values = c("Fixed" = 15,"Varied" = 17)) +
  # labs(x=TeX("$\\gamma$"), y="Bias (%)", title = "(III)") +
  # labs(x=TeX("$\\gamma$"), y="Abs. bias", title = "(III)") +
  labs(x=TeX("$\\gamma$"), y="Abs. bias") +
  guides(fill = guide_legend(override.aes = list(size = 12))) +
  theme_pubclean() +
  theme(axis.text.x = element_text(size =26, face = "bold"),
        axis.text.y = element_text(size =26, face = "bold"),
        # axis.text.y = element_blank(),
        axis.title.x = element_text(size = 44, face = "bold"),
        axis.title.y = element_text(size = 26, face="bold"),
        # axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=26, face = "bold"),
        legend.position = "top",
        # legend.position = c(0.1,0.9),
        # legend.position = "none",
        # legend.direction = "horziontal",
        legend.title = element_blank(),
        legend.text = element_text(size = 28, face = "bold"),
        # legend.key.size = unit(1.2,"cm"),
        plot.title = element_text(size = 36, face = "bold", hjust = .5, vjust = -0.5))

legend.plot <- ggplot(data.frame(name=c("HT","Hajek","Fixed","Varied"),
                        x=seq(4),y=seq(4)),
             aes(x = x, y = x,
                 color = factor(name,levels = c("HT","Hajek","Fixed","Varied")),
                 shape = factor(name,levels = c("HT","Hajek","Fixed","Varied"))))+
  geom_point()+
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  scale_color_manual(values = c("HT" = "#990000","Hajek" = "#0065A9",
                                "Fixed" = "#009933", "Varied" = "#660033"),
                     labels = c("HT","Hajek","Fixed","Varied")) +
  scale_shape_manual(values = c("HT" = 16,"Hajek" = 18,
                                "Fixed" = 15,"Varied" = 17),
                     labels = c("HT","Hajek","Fixed","Varied")) +
  guides(colour = guide_legend(override.aes = list(size=12))) +
  theme(legend.title = element_blank(),
        legend.position = c(0.5,0.5),
        legend.text = element_text(size = 28, face = "bold"),
        legend.key.size = unit(1.2,"cm"),
        legend.spacing = unit(5,"cm"))

plot.contamination.legend <- ggarrange(bias_cluster_plot, legend.plot,
                              nrow = 1, ncol = 2, widths = c(1.75,1), align = "h") +
                                theme(plot.margin=unit(c(-0.1,1,1,1),"cm"))

plot.all.combined <- ggarrange(plot.noise.censored.combined, plot.contamination.legend,
                               nrow = 2)
# 
# ggsave(filename = "Reproducibility/Simulations/graphics/Main/Bias_plot_all_scenarios_combined_1000iter.jpeg",
#        plot = plot.all.combined,
#        width = 20, height = 12)
ggsave(filename = "Reproducibility/Simulations/graphics/Appendix/Bias_cross_clusters_contamination.jpeg",
       plot = bias_cluster_plot,
       width = 16, height = 8)



