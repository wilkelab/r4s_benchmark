library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(readr)

setwd("r4s_benchmark/")
r_nobias_dNdS_raw <- read_csv("mech_codon/processed_rates/all_r4s_orig_rates_gamma_nobias.csv")

#####################################################################
### Calculating correlations, RMSD, and bias                      ###
#####################################################################

gamma_t <- data.frame(gamma_distr=c(1:6),
                      shape=c(0.23, 0.35, 0.647, 0.312, 0.378, 0.238),
                      rate=c(0.029, 0.024, 0.070, 0.032, 0.062, 0.018))

r_nobias_dNdS <- r_nobias_dNdS_raw %>%
  left_join(gamma_t) %>%
  na.omit() %>%
  group_by(bl,num_taxa,rep,gamma_distr) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true)) %>% 
  summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
            rmsd_true=sqrt(mean((score_norm - true_norm)^2)))

#################################################################################
### Plot 1.                                                                   ###
### Correlations between Rate4Site vs true and inferred dN/dS for simulations ###
### without codon bias and gamma distributed true rates.                     ###
### RMSD between Rate4Site vs dN/dS no codon bias                             ###
#################################################################################

plot_lst <- list()
for (i in c(1:6)) {
  r <- r_nobias_dNdS %>% filter(gamma_distr==i)
  
  ############ Correlations vs Branch len ###########
  colfunc <- colorRampPalette(c("cyan2","navyblue"))
  cor_p <- ggplot(r,aes(bl,cor_true,colour=factor(num_taxa)))+
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.4)+
    scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
    scale_colour_manual(values=colfunc(5)) +
    guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
    stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
    xlab("Branch Length") +
    ylab("Correlation (spearman)") +
    coord_cartesian(ylim=c(0, 1),xlim=c(0.0023,0.66))+
    scale_y_continuous(breaks=seq(0,1,0.2))+
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12),
          legend.position="none")
  
  ############### RMSD vs Branch len ###############
  rmsd_p <- ggplot(r,aes(bl,rmsd_true,colour=factor(num_taxa))) + 
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
                 fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
                 geom = "pointrange",
                 size=0.4)+
    scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
    scale_colour_manual(values=colfunc(5)) +
    guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
    stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
    xlab("Branch Length") +
    ylab("RMSD") +
    coord_cartesian(ylim=c(0,2.5),xlim=c(0.0023,0.66))+
    scale_y_continuous(breaks=seq(0,2.5,0.5))+
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 12),
          legend.position="none")
  
  plot_lst[[length(plot_lst)+1]] <- cor_p
  plot_lst[[length(plot_lst)+1]] <- rmsd_p
}
################### Plotting the grid ###################
grobs <- ggplotGrob(cor_p)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow <- plot_grid(plotlist=plot_lst,
                  labels="AUTO",
                  align = 'vh',
                  hjust = -1,
                  ncol=2,
                  nrow=6)

p <- plot_grid(prow, legend, rel_widths = c(2, .3))

save_plot("plots/r4s_v_dNdS_gamma.png", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 6, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)
