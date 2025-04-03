###
# Script that generates the graphic results of the NMR bias-variance tradeoff simulations
###


# Load --------------------------------------------------------------------


library(ggplot2)
library(data.table)
library(latex2exp)
library(ggpubr)



# Joint -------------------------------------------------------------------


# csv_names <- dir("Reproducibility/Simulations/results/NMR_PA/")
# path_name <- "Reproducibility/Simulations/results/NMR_PA/"
# 
# combined_dt <- data.table()
# 
# for (res_name in csv_names) {
# 
#   curr_dt <- fread(paste0(path_name,res_name))
#   combined_dt <- rbindlist(list(combined_dt,curr_dt))
# }
# 
# 
# setorderv(combined_dt, c("K","adj.mat.used","iter"))
# 
# write.csv(combined_dt,
#           "Reproducibility/Simulations/results/MR_bias_var_PA.csv",
#           row.names = FALSE)
# 



# Noisy networks ----------------------------------------------------------------

MR.sim.results <- fread("Reproducibility/Simulations/results/MR_bias_var_PA.csv")

MR.sim.results[,true_effect := rep(c(1,0.75,0.5,0.25,0.5),nrow(MR.sim.results)/5)]

MR.sim.results[,hajek_coverage := ifelse(true_effect >= hajek_ce - 1.96*sqrt(var_hajek_ce) &
                              true_effect <= hajek_ce + 1.96*sqrt(var_hajek_ce),1,0)]

MR.sim.results[,ht_coverage := ifelse(true_effect >= ht_ce - 1.96*sqrt(var_ht_ce) &
                              true_effect <= ht_ce + 1.96*sqrt(var_ht_ce),1,0)]

mean.sim.results <- MR.sim.results[,.(m_ht = abs(mean(ht_ce-true_effect)),
                                   m_hajek = abs(mean(hajek_ce-true_effect)),
                                   se_ht = sd(ht_ce,na.rm=T),
                                   se_hajek = sd(hajek_ce,na.rm=T),
                                   rmse_ht = sqrt((mean(ht_ce-true_effect))^2 +
                                                    var(ht_ce,na.rm=T)),
                                   rmse_hajek = sqrt((mean(hajek_ce-true_effect))^2 +
                                                       var(hajek_ce,na.rm=T)),
                                   coverage_ht = mean(ht_coverage, na.rm=T),
                                   coverage_hajek = mean(hajek_coverage, na.rm=T)
                                   ),
                                by = c("ce_contrast","K","adj.mat.used","with.true.adj")]



# Main plot c11-c10

p.hajek.bias <- ggplot(
                    mean.sim.results[ce_contrast == "c11-c10",],
                         aes(x=factor(K),y = m_hajek,
                             col = with.true.adj, shape = with.true.adj)) +
                   labs(title = "Bias", 
                         x = "",
                         y = ""
                         ) +
                  geom_point(size = 10, stroke = 2) +
                  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
                  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
                  scale_y_continuous(breaks = seq(0,0.4,0.1),limits = c(0,0.4)) +
                  theme_pubclean() +
                  guides(color = guide_legend(override.aes = list(size = 12))) +
                 theme(axis.title.x = element_text(size = 28, face = "bold"),
                  title = element_text(size = 26, face="bold"),
                  axis.text = element_text(size=26, face = "bold"),
                  legend.title = element_blank(),
                  legend.text = element_text(face = "bold", size = 28),
                  legend.key.size = unit(1.2,"cm"),
                  plot.margin=unit(c(1,0.15,1,1),"cm")) 


                
              
p.hajek.sd <- ggplot(mean.sim.results[ce_contrast == "c11-c10",],
                     aes(x=factor(K),y = se_hajek,
                         col = with.true.adj, shape = with.true.adj)) +
                labs(title = "SD", 
                     # x = "",
                     x = "# networks used",
                     # y = TeX('$\\tau - \\hat{\\tau}$')
                     y = ""
                     ) +
              geom_point(size = 10, stroke = 2) +
              scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
              scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
              scale_y_continuous(breaks = seq(0,0.4,0.1),limits = c(0,0.4)) +
              theme_pubclean() +
              guides(color = guide_legend(override.aes = list(size = 12))) +
              theme(axis.title.x = element_text(size = 28, face = "bold"),
                    title = element_text(size = 26, face="bold"),
                    axis.text = element_text(size=26, face = "bold"),
                    axis.text.y = element_blank(),
                    legend.title = element_blank(),
                    legend.text = element_text(face = "bold", size = 28),
                    legend.key.size = unit(1.2,"cm"),
                    plot.margin=unit(c(1,0.15,1,0.15),"cm")) 


p.hajek.rmse <- ggplot(mean.sim.results[ce_contrast == "c11-c10",],
                     aes(x=factor(K),y = rmse_hajek,
                         col = with.true.adj, shape = with.true.adj)) +
                labs(title = "RMSE", 
                     # x = "# networks used",
                     x = "",
                     # y = TeX('$\\tau - \\hat{\\tau}$')
                     y = ""
                     ) +
                geom_point(size = 10, stroke = 2) +
                scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
                scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
                scale_y_continuous(breaks = seq(0,0.4,0.1),limits = c(0,0.4)) +
                theme_pubclean() +
                guides(color = guide_legend(override.aes = list(size = 12))) +
                theme(axis.title.x = element_text(size = 28, face = "bold"),
                      title = element_text(size = 26, face="bold"),
                      axis.text = element_text(size=26, face = "bold"),
                      # legend.title = element_text(size = 14),
                      axis.text.y = element_blank(),
                      legend.title = element_blank(),
                      legend.text = element_text(face = "bold", size = 28),
                      legend.key.size = unit(1.2,"cm"),
                      plot.margin=unit(c(1,1,1,0.15),"cm")) 
  

p.hajek.both <- ggarrange(p.hajek.bias,p.hajek.sd,p.hajek.rmse,
                          nrow = 1,ncol = 3,
                          widths = c(1,1,1),
                          common.legend = T, legend = "top", align = "h") 

ggsave("Reproducibility/Simulations/graphics/Main/MR_bias_var_PA_n3000_eta025_c11_c10.jpeg",
       p.hajek.both,
       dpi = 700,
       width = 20, height = 10)



p.hajek.coverage <- ggplot(mean.sim.results[ce_contrast == "c11-c10",],
                       # aes(x=factor(K),y = coverage_hajek,
                       aes(x=factor(K),y = coverage_ht,
                           col = with.true.adj, shape = with.true.adj)) +
  labs(title = "coverage", 
       # x = "# networks used",
       x = "",
       # y = TeX('$\\tau - \\hat{\\tau}$')
       y = ""
  ) +
  geom_point(size = 10, stroke = 2) +
  geom_line() + 
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  geom_abline(intercept = 0.95, slope = 0, lty = "dashed", linewidth = 1.1) +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1)) +
  theme_pubclean() +
  guides(color = guide_legend(override.aes = list(size = 12))) +
  theme(axis.title.x = element_text(size = 28, face = "bold"),
        title = element_text(size = 26, face="bold"),
        axis.text = element_text(size=26, face = "bold"),
        # legend.title = element_text(size = 14),
        # axis.text.y = element_blank(),
        # legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 28),
        legend.key.size = unit(1.2,"cm"),
        plot.margin=unit(c(1,1,1,0.15),"cm")) 



mean.sim.results$labels = factor(mean.sim.results$ce_contrast,
                                    labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))

# APPENDIX plot

p.hajek.bias.apdx <- 
  ggplot(mean.sim.results[ce_contrast %in% c("c01-c00","c11-c00"),],
  # ggplot(mean.sim.results,
                       aes(x=factor(K),y = m_hajek,
                           col = with.true.adj, shape = with.true.adj)) +
  labs(title = "Bias", 
       x = "",
       y = ""
  ) +
  geom_point(size = 8, stroke = 1.6) +
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  scale_y_continuous(breaks = seq(0,0.3,0.1),limits = c(0,0.3)) +
  facet_wrap(~labels, nrow = 1, labeller = label_parsed) +
  theme_pubclean() +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"),
        title = element_text(size = 26, face="bold"),
        axis.text = element_text(size=22, face = "bold"),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        legend.key.size = unit(1.1,"cm"),
        strip.text = element_text(size=26, face = "bold")) 



p.hajek.sd.apdx <- 
  ggplot(mean.sim.results[ce_contrast %in% c("c01-c00","c11-c00"),],
  # ggplot(mean.sim.results,
                     aes(x=factor(K),y = se_hajek,
                         col = with.true.adj, shape = with.true.adj)) +
  labs(title = "SD", 
       x = "",
       y = ""
  ) +
  geom_point(size = 8, stroke = 1.6) +
 
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  scale_y_continuous(breaks = seq(0,0.3,0.1),limits = c(0,0.3)) +
  facet_wrap(~labels, nrow = 1, labeller = label_parsed) +
  theme_pubclean() +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"),
        title = element_text(size = 26, face="bold"),
        axis.text = element_text(size=22, face = "bold"),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        legend.key.size = unit(1.1,"cm"),
        strip.text = element_text(size=26, face = "bold")) 



p.hajek.rmse.apdx <- 
  ggplot(mean.sim.results[ce_contrast %in% c("c01-c00","c11-c00"),],
                       aes(x=factor(K),y = rmse_hajek,
                           col = with.true.adj, shape = with.true.adj)) +
                    labs(title = "RMSE", 
                         x = "# networks used",
                         y = ""
                    ) +
                    geom_point(size = 8, stroke = 1.6) +
                    scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
                    scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
                    scale_y_continuous(breaks = seq(0,0.3,0.1),limits = c(0,0.3)) +
                    facet_wrap(~labels, nrow = 1, labeller = label_parsed) +
                    theme_pubclean() +
                    guides(color = guide_legend(override.aes = list(size = 10))) +
                    theme(axis.title.x = element_text(size = 26, face = "bold"),
                          title = element_text(size = 26, face="bold"),
                          axis.text = element_text(size=22, face = "bold"),
                          legend.title = element_blank(),
                          legend.text = element_text(face = "bold", size = 26),
                          legend.key.size = unit(1.1,"cm"),
                          strip.text = element_text(size=26, face = "bold")) 
                  
                  



p.hajek.both.apdx <- ggarrange(p.hajek.bias.apdx,p.hajek.sd.apdx,p.hajek.rmse.apdx,
                          nrow = 3,ncol = 1,
                          # widths = c(1,1,2),
                          common.legend = T, legend = "top")

ggsave("Reproducibility/Simulations/graphics/Appendix/MR_bias_var_PA_n3000_eta025_APDX.jpeg",
       p.hajek.both.apdx,
       height = 24,
       width = 18)


# SE/SD comparisons -------------------------------------------------------


se_sd_comparison <- MR.sim.results[with.true.adj==TRUE ,
                                   .(hajek.se.sd = sqrt(var_hajek_ce)/sd(hajek_ce),
                                       ht.se.sd = sqrt(var_ht_ce)/sd(ht_ce)),
                                   by = c("adj.mat.used","ce_contrast","K")]

se_sd_comparison <- melt.data.table(se_sd_comparison,
                      id.vars = c("adj.mat.used","ce_contrast","K"),
                      measure.vars = c("hajek.se.sd","ht.se.sd"),
                      variable.name = "esti",
                      value.name = "se.sd")

se_sd_comparison_mean <- MR.sim.results[with.true.adj==TRUE ,
                                   .(hajek.se.sd = mean(sqrt(var_hajek_ce)/sd(hajek_ce),na.rm=T),
                                       ht.se.sd = mean(sqrt(var_ht_ce)/sd(ht_ce),na.rm=T)),
                                   by = c("adj.mat.used","ce_contrast","K")]

se_sd_comparison_mean <- melt.data.table(se_sd_comparison_mean,
                      id.vars = c("adj.mat.used","ce_contrast","K"),
                      measure.vars = c("hajek.se.sd","ht.se.sd"),
                      variable.name = "esti",
                      value.name = "se.sd")

# [ce_contrast=="c01-c00",]
se.to.sd.plot <- ggplot(
  # se_sd_comparison_mean[ce_contrast=="c01-c00",],
  se_sd_comparison_mean[ce_contrast=="c11-c00",],
               aes(x=factor(K), y = se.sd, col=esti, shape = esti)) +
            geom_point(size=10,alpha= 0.7) +
            geom_hline(yintercept = 1, lty = "dashed", linewidth = 1.1) +
            scale_color_manual(values = c("ht.se.sd" = "#990000","hajek.se.sd" = "#0065A9")
                               ,labels = c("Hajek","HT")
                               ) +
            scale_shape_manual(values = c("ht.se.sd" = 16,"hajek.se.sd" = 18)
                               ,labels = c("Hajek","HT")
                               ) +
            labs(x="# networks used",
                 y = "SE/SD") +
          guides(col = guide_legend(override.aes = list(size = 12))) +
          theme_pubclean() +
          theme(axis.text.x = element_text(size =20, face = "bold"),
                axis.text.y = element_text(size =20, face = "bold"),
                axis.title.x = element_text(size = 26, face="bold", vjust = 0.2),
                axis.title.y = element_text(size = 22, face="bold"),
                legend.position = "top",
                legend.title = element_blank(),
                legend.text = element_text(size = 26, face = "bold"),
                legend.key.size = unit(1.2,"cm"),
          )

ggsave("Reproducibility/Simulations/graphics/Appendix/MR_PA_n3000_eta025_c01_c00_SE_to_SD_APDX.jpeg",
       se.to.sd.plot,
       height = 10,
       width = 14)

# Palluck et al. Networks -------------------------------------------------

MR.sim.palluck <- fread("Reproducibility/Simulations/results/MR_bias_var_PalluckNets.csv")

MR.sim.palluck[,true_effect := rep(c(1,0.75,0.5,0.25,0.5),nrow(MR.sim.palluck)/5)]

mean.sim.palluck <- MR.sim.palluck[,.(m_ht = abs(mean(ht_ce-true_effect)),
                                      m_hajek = abs(mean(hajek_ce-true_effect)),
                                      se_ht = sd(ht_ce,na.rm=T),
                                      se_hajek = sd(hajek_ce,na.rm=T),
                                      rmse_ht = sqrt((mean(ht_ce-true_effect))^2 +
                                                       var(ht_ce,na.rm=T)),
                                      rmse_hajek = sqrt((mean(hajek_ce-true_effect))^2 +
                                                          var(hajek_ce,na.rm=T))),
                                   by = c("ce_contrast","K","adj.mat.used","with.true.adj")]

mean.sim.palluck$labels = factor(mean.sim.palluck$ce_contrast,
                                    labels = c(TeX(r"($\tau(c_{\0\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\0}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\0})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\0\1})$)",output = "character"),
                                               TeX(r"($\tau(c_{\1\1}, c_{\1\0})$)",output = "character")))


p.hajek.bias.apdx.palluck <- 
  ggplot(mean.sim.palluck[ce_contrast %in% c("c11-c00","c11-c10"),],
         aes(x=factor(K),y = m_hajek,
             col = with.true.adj, shape = with.true.adj)) +
  labs(title = "Bias", 
       x = "",
       # y = TeX('$\\tau - \\hat{\\tau}$')
       y = ""
  ) +
  geom_point(size = 8, stroke = 1.6) +
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  # scale_y_continuous(expand = expansion(add = .15)) +
  scale_y_continuous(breaks = seq(0,0.6,0.1),limits = c(0,0.6)) +
  facet_wrap(~labels,nrow = 1, labeller = label_parsed, scales = "free_y") +
  theme_pubclean() +
  # theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"),
        title = element_text(size = 26, face="bold"),
        axis.text = element_text(size=22, face = "bold"),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        legend.key.size = unit(1.1,"cm"),
        strip.text = element_text(size=22, face = "bold")) 



p.hajek.sd.apdx.palluck <- 
  ggplot(mean.sim.palluck[ce_contrast %in% c("c11-c00","c11-c10"),],
         aes(x=factor(K),y = se_hajek,
             col = with.true.adj, shape = with.true.adj)) +
  labs(title = "SD", 
       x = "",
       # y = TeX('$\\tau - \\hat{\\tau}$')
       y = ""
  ) +
  geom_point(size = 8, stroke = 1.6) +
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  # scale_y_continuous(expand = expansion(add = .15)) +
  scale_y_continuous(breaks = seq(0,0.6,0.1),limits = c(0,0.6)) +  facet_wrap(~labels,nrow = 1, labeller = label_parsed, scales = "free_y") +
  theme_pubclean() +
  # theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"),
        title = element_text(size = 26, face="bold"),
        axis.text = element_text(size=22, face = "bold"),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        legend.key.size = unit(1.1,"cm"),
        strip.text = element_text(size=22, face = "bold")) 



p.hajek.rmse.apdx.palluck <- 
  ggplot(mean.sim.palluck[ce_contrast %in% c("c11-c00","c11-c10"),],
         aes(x=factor(K),y = rmse_hajek,
             col = with.true.adj, shape = with.true.adj)) +
  labs(title = "RMSE", 
       x = "# networks used",
       # y = TeX('$\\tau - \\hat{\\tau}$')
       y = ""
  ) +
  geom_point(size = 8, stroke = 1.6) +
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  # scale_y_continuous(expand = expansion(add = .15)) +
  scale_y_continuous(breaks = seq(0,0.6,0.1),limits = c(0,0.6)) +  facet_wrap(~labels,nrow = 1, labeller = label_parsed, scales = "free_y") +
  theme_pubclean() +
  # theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 10))) +
  theme(axis.title.x = element_text(size = 26, face = "bold"),
        title = element_text(size = 26, face="bold"),
        axis.text = element_text(size=22, face = "bold"),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 26),
        legend.key.size = unit(1.1,"cm"),
        strip.text = element_text(size=22, face = "bold")) 



p.hajek.all.apdx.palluck <- ggarrange(p.hajek.bias.apdx.palluck,
                               p.hajek.sd.apdx.palluck,
                               p.hajek.rmse.apdx.palluck,
                               nrow = 3,ncol = 1,
                               # widths = c(1,1,2),
                               common.legend = T, legend = "top")

ggsave("Reproducibility/Simulations/graphics/Appendix/NMR_bias_var_APDX_Palluck_Networks.jpeg",
       plot = p.hajek.all.apdx.palluck,
       width = 18, height = 24)




