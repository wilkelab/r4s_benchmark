library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

model="mech_codon"
r_bias <- read.csv(paste0(model,"/r4s_rates/processed_rates/all_r4s_norm_rates_bias.csv"),row.names=1)
r_nobias <- read.csv(paste0(model,"/r4s_rates/processed_rates/all_r4s_norm_rates_nobias.csv"),row.names=1)

##############################################################
### Plotting correlations between true rates and r4s rates ###
##############################################################

####################### Bias #################################
colfunc <- colorRampPalette(c("orange","darkred"))
p_bias_num_taxa_cor_true <- ggplot(r_bias,aes(num_taxa,cor_true,color=factor(bl),group=bl)) + 
  stat_summary(fun.y = mean,
              fun.ymin = function(x) mean(x) - sd(x), 
             fun.ymax = function(x) mean(x) + sd(x), 
              geom = "pointrange",
             size=0.6) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Branch Length",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",size=0.6,aes(color=factor(bl)))+
  xlab("Number of Taxa") +
  ylab("Correlation (spearman)") +
	scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
  scale_x_log10(breaks=c(128,256,512,1024,2048)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_bias_cor_true_v_num_taxa.png"))

colfunc <- colorRampPalette(c("cyan2","navyblue"))
p_bias_bl_cor_true <- ggplot(r_bias,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.6)+
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
  xlab("Branch Length") +
  ylab("Correlation (spearman)") +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_bias_cor_true_v_branch_len.png"))

####################### No bias #################################

colfunc <- colorRampPalette(c("orange","darkred"))
p_nobias_num_taxa_cor_true <- ggplot(r_nobias,aes(num_taxa,cor_true,color=factor(bl),group=bl)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.6) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Branch Length",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",size=0.6,aes(color=factor(bl)))+
  xlab("Number of Taxa") +
  ylab("Correlation (spearman)") +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
  scale_x_log10(breaks=c(128,256,512,1024,2048)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_nobias_cor_true_v_num_taxa.png"))

colfunc <- colorRampPalette(c("cyan2","navyblue"))
p_nobias_bl_cor_true <- ggplot(r_nobias,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.6) +
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",size=0.6,aes(color=factor(num_taxa)))+
  xlab("Branch Length") +
  ylab("Correlation (spearman)") +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_nobias_cor_true_v_branch_len.png"))

##################################################################
### Plotting correlations between assigned rates and r4s rates ###
##################################################################

####################### Bias #################################
colfunc <- colorRampPalette(c("orange","darkred"))
p_bias_num_taxa_cor_assigned <- ggplot(r_bias,aes(num_taxa,cor_assigned,color=factor(bl),group=bl)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.6) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Branch Length",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",size=0.6,aes(color=factor(bl)))+
  xlab("Number of Taxa") +
  ylab("Correlation (spearman)") +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
  scale_x_log10(breaks=c(128,256,512,1024,2048)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_bias_cor_assigned_v_num_taxa.png"))

colfunc <- colorRampPalette(c("cyan2","navyblue"))
p_bias_bl_cor_assigned <- ggplot(r_bias,aes(bl,cor_assigned,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.6)+
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
  xlab("Branch Length") +
  ylab("Correlation (spearman)") +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_bias_cor_assigned_v_branch_len.png"))

#######################No bias #################################
colfunc <- colorRampPalette(c("orange","darkred"))
p_nobias_num_taxa_cor_assigned <- ggplot(r_nobias,aes(num_taxa,cor_assigned,color=factor(bl),group=bl)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.6) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Branch Length",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",size=0.6,aes(color=factor(bl)))+
  xlab("Number of Taxa") +
  ylab("Correlation (spearman)") +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
  scale_x_log10(breaks=c(128,256,512,1024,2048)) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_nobias_cor_assigned_v_num_taxa.png"))

colfunc <- colorRampPalette(c("cyan2","navyblue"))
p_nobias_bl_cor_assigned <- ggplot(r_nobias,aes(bl,cor_assigned,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.6) +
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",size=0.6,aes(color=factor(num_taxa)))+
  xlab("Branch Length") +
  ylab("Correlation (spearman)") +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
ggsave(paste0("plots/",model,"_nobias_cor_assigned_v_branch_len.png"))

##############################################################
### Plotting combined plots                                ###
##############################################################

pg_bias <- plot_grid(p_bias_num_taxa_cor_true,
                 p_bias_num_taxa_cor_assigned,
                 p_bias_bl_cor_true,
                 p_bias_bl_cor_assigned,
                 labels = c("True rate (bias)", "Assigned rate (bias)"),
                 ncol=2,nrow=2)
save_plot(paste0("plots/",model,"_bias_combined.png"), pg,
           ncol = 2, # we're saving a grid plot of 2 columns
           nrow = 2, # and 2 rows
           # each individual subplot should have an aspect ratio of 1.3
           base_aspect_ratio = 1.3)

pg_nobias <- plot_grid(p_nobias_num_taxa_cor_true,
                p_nobias_num_taxa_cor_assigned,
                p_nobias_bl_cor_true,
                p_nobias_bl_cor_assigned,
                labels = c("True rate (no bias)", "Assigned rate (no bias)"),
                ncol=2,nrow=2)
save_plot(paste0("plots/",model,"_nobias_combined.png"), pg,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)


