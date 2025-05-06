###
# Script that plot the Number of Effective Units (NEU) of the NMR estimator in the simulations
###

# Load --------------------------------------------------------------------


library(ggplot2)
library(latex2exp)
library(ggpubr)
library(tidyverse)  
library(data.table)

# Merge files ---------------------------------------------------------------

# csv_names <- dir("Reproducibility/Simulations/neu_paluck/")
# path_name <- "Reproducibility/Simulations/neu_paluck/"
# csv_names <- dir("Reproducibility/Simulations/neu/")
# path_name <- "Reproducibility/Simulations/neu/"

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
#           "Reproducibility/Simulations/results/NMR_paluck_neu_results.csv",
#           row.names = FALSE)
# write.csv(combined_dt,
#           "Reproducibility/Simulations/results/NMR_neu_results.csv",
#           row.names = FALSE)



# Load the data ------------------------------------------------------------

neu_raw <- read_csv("Reproducibility/Simulations/results/NMR_neu_results.csv")


# ---- Summarise across iterations ----
neu_summary <- neu_raw %>%
  group_by(K, with.true.adj, exposure) %>%
  summarise(mean_neu = mean(NEU),
            sd_neu   = sd(NEU),
            .groups  = "drop")

# ---- Plot ----
pd <- position_dodge(width = 0.3)     

neu_nmr_plot <- ggplot(neu_summary,
       aes(x = factor(K),
           y = mean_neu,
           color = with.true.adj,
           shape  = with.true.adj)) +
  geom_point(size = 3, position = pd, stroke=1.5) +
  geom_errorbar(aes(ymin = mean_neu - sd_neu,
                    ymax = mean_neu + sd_neu),
                width = 0.15,
                linewidth = .7,
                position = pd) +
  # scale_colour_manual(values = c(`TRUE`  = "green4",
  #                                `FALSE` = "red4"),
  #                     name = "Includes true network") +
  # scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 1),
  #                    name = "Includes true network") +
  scale_color_manual(labels = c("A* not included", "A* included"), values = c("red4","green4")) +
  scale_shape_manual(labels = c("A* not included", "A* included"), values = c(1,4)) +
  facet_wrap(~ exposure, nrow = 2,
             labeller = labeller(exposure = function(x) paste("Exposure =", x)),
             scales = "fixed") +           
  labs(x = "# networks used",
       y = "Mean NEU (± 1 SD)",
       title = "",) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text = element_text(size=16, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold", size = 16),
        legend.key.size = unit(1.2,"cm"),
        legend.position = "top",)


ggsave("Reproducibility/Simulations/graphics/Appendix/NMR_NEU_plot.jpeg",
       plot = neu_nmr_plot,
       width = 14, height = 8)

