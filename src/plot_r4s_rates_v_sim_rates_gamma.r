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

r_nobias_dNdS <- r_nobias_dNdS_raw %>%
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>% 
  summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
            rmsd_true=sqrt(mean((score_norm - true_norm)^2)),
            cor_inferred=cor(score,inferred,method="spearman",use="pairwise.complete.obs"),
            rmsd_inferred=sqrt(mean((score_norm - inferred_norm)^2)))

#################################################################################
### Plot 1.                                                                   ###
### Correlations between Rate4Site vs true and inferred dN/dS for simulations ###
### without codon bias and gamma distributed true rates.                     ###
### RMSD between Rate4Site vs dN/dS no codon bias                             ###
#################################################################################

############ Correlations vs Branch len ###########
colfunc <- colorRampPalette(c("cyan2","navyblue"))
dNdS_p_nobias_bl_cor_true <- ggplot(r_nobias_dNdS,aes(bl,cor_true,colour=factor(num_taxa)))+
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
        legend.title = element_text(size = 12))

############### RMSD vs Branch len ###############
dNdS_p_nobias_bl_rmsd_true <- ggplot(r_nobias_dNdS,aes(bl,rmsd_true,colour=factor(num_taxa))) + 
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
  coord_cartesian(ylim=c(0,5),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,5,1))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

############ Correlations vs Branch len ###########
colfunc <- colorRampPalette(c("cyan2","navyblue"))
dNdS_p_nobias_bl_cor_inferred <- ggplot(r_nobias_dNdS,aes(bl,cor_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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
        legend.title = element_text(size = 12))

############### RMSD vs Branch len ###############
dNdS_p_nobias_bl_rmsd_inferred <- ggplot(r_nobias_dNdS,aes(bl,rmsd_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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
  coord_cartesian(ylim=c(0,5),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,5,1))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(dNdS_p_nobias_bl_cor_true)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_bias <- plot_grid(dNdS_p_nobias_bl_cor_true+theme(legend.position="none")+ggtitle("True dN/dS"),
                       dNdS_p_nobias_bl_cor_inferred+theme(legend.position="none",axis.title.y = element_blank())+ggtitle("Inferred dN/dS"),
                       dNdS_p_nobias_bl_rmsd_true+theme(legend.position="none"),
                       dNdS_p_nobias_bl_rmsd_inferred+theme(legend.position="none",axis.title.y = element_blank()),
                       labels=c("A","C","B","D"),
                       align = 'vh',
                       hjust = -1,
                       ncol=2,
                       nrow=2)

p <- plot_grid( prow_bias, legend, rel_widths = c(2, .3))

save_plot("plots/r4s_v_dNdS_true_gamma.png", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)
