library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)

setwd("r4s_benchmark/")
r_bias_dNdS_raw <- read.csv("mech_codon/r4s_rates/processed_rates/all_r4s_orig_rates_bias.csv",row.names=1)
r_nobias_dNdS_raw <- read.csv("mech_codon/r4s_rates/processed_rates/all_r4s_orig_rates_nobias.csv",row.names=1)
r_bias_MutSel_raw <- read.csv("mut_sel/r4s_rates/processed_rates/all_r4s_orig_rates_bias.csv",row.names=1)
r_nobias_MutSel_raw <- read.csv("mut_sel/r4s_rates/processed_rates/all_r4s_orig_rates_nobias.csv",row.names=1)

#####################################################################
### Calculating correlations, RMSD, and bias                      ###
#####################################################################

sim1 <- r_bias_dNdS_raw %>% group_by(bl,num_taxa,rep)
sim1 %>% summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
                          rmsd_true=sqrt(sum((score_norm - true_norm)^2)),
                          bias_true=sum(score_norm - true_norm)) -> r_bias_dNdS

sim2 <- r_nobias_dNdS_raw %>% group_by(bl,num_taxa,rep)
sim2 %>% summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
                          rmsd_true=sqrt(sum((score_norm - true_norm)^2)),
                          bias_true=sum(score_norm - true_norm)) -> r_nobias_dNdS


sim3 <- r_bias_MutSel_raw %>% group_by(bl,num_taxa,rep)
sim3 %>% summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
                   rmsd_true=sqrt(sum((score_norm[!is.na(true_norm)] - true_norm[!is.na(true_norm)])^2)),
                   bias_true=sum(score_norm[!is.na(true_norm)] - true_norm[!is.na(true_norm)]),
                   cor_inferred=cor(score,inferred,method="spearman",use="pairwise.complete.obs"),
                   rmsd_inferred=sqrt(sum((score_norm[!is.na(inferred_norm)] - inferred_norm[!is.na(inferred_norm)])^2)),
                   bias_inferred=sum(score_norm[!is.na(inferred_norm)] - inferred_norm[!is.na(inferred_norm)])) -> r_bias_MutSel

sim4 <- r_nobias_MutSel_raw %>% group_by(bl,num_taxa,rep)
sim4 %>% summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
                   rmsd_true=sqrt(sum((score_norm[!is.na(true_norm)] - true_norm[!is.na(true_norm)])^2)),
                   bias_true=sum(score_norm[!is.na(true_norm)] - true_norm[!is.na(true_norm)]),
                   cor_inferred=cor(score,inferred,method="spearman",use="pairwise.complete.obs"),
                   rmsd_inferred=sqrt(sum((score_norm[!is.na(inferred_norm)] - inferred_norm[!is.na(inferred_norm)])^2)),
                   bias_inferred=sum(score_norm[!is.na(inferred_norm)] - inferred_norm[!is.na(inferred_norm)])) -> r_nobias_MutSel

#####################################################################
### Plotting Figure 2:                                            ###
### correlations between Rate4Site vs dN/dS/MutSel no codon bias  ###
### RMSD between Rate4Site vs dN/dS/MutSel no codon bias          ###
### bias between Rate4Site vs dN/dS/MutSel no codon bias          ###
#####################################################################

################### Correlations ###################
colfunc <- colorRampPalette(c("cyan2","navyblue"))
dNdS_p_nobias_bl_cor_true <- ggplot(r_nobias_dNdS,aes(bl,cor_true,colour=factor(num_taxa)))+
  ggtitle("dN/dS") +
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
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_nobias_bl_cor_true <- ggplot(r_nobias_MutSel,aes(bl,cor_true,colour=factor(num_taxa)))+ 
  ggtitle("MutSel") +
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
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))  

################### RMSD ###################
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
  scale_y_continuous(breaks=seq(0,16,2), limits = c(0,16))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_nobias_bl_rmsd_true <- ggplot(r_nobias_MutSel,aes(bl,rmsd_true,colour=factor(num_taxa)))+ 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
               geom = "pointrange",
               size=0.4,pseudocount=0.0001)+
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
  xlab("Branch Length") +
  ylab("RMSD") +
  scale_y_continuous(breaks=seq(0,16,2), limits = c(0,16))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### bias ###################
dNdS_p_nobias_bl_bias_true <- ggplot(r_nobias_dNdS,aes(bl,bias_true,colour=factor(num_taxa))) + 
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
  ylab("Bias") +
  scale_y_continuous(breaks=c(-3:0),limits = c(-3,0.25))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_nobias_bl_bias_true <- ggplot(r_nobias_MutSel,aes(bl,bias_true,colour=factor(num_taxa)))+
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
  ylab("Bias") +
  scale_y_continuous(breaks=c(-3:0),limits = c(-4,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(MutSel_p_nobias_bl_bias_true)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_nobias <- plot_grid(dNdS_p_nobias_bl_cor_true+theme(legend.position="none",axis.title.x = element_blank()),
                MutSel_p_nobias_bl_cor_true+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
                dNdS_p_nobias_bl_rmsd_true+theme(legend.position="none",axis.title.x = element_blank()),
                MutSel_p_nobias_bl_rmsd_true+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
                dNdS_p_nobias_bl_bias_true+theme(legend.position="none"),
                MutSel_p_nobias_bl_bias_true+theme(legend.position="none",axis.title.y = element_blank()),
                align = 'vh',
                hjust = -1,
                ncol=2,
                nrow=3)
p_nobias <- plot_grid( prow_nobias, legend, rel_widths = c(2, .3))
save_plot("plots/fig2_r4s_v_dNdS_MutSel_nobias.png", p_nobias,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

######################################################################
### Plotting Figure 3:                                             ###
### correlations between Rate4Site vs dN/dS/MutSel with codon bias ###
### RMSD between Rate4Site vs dN/dS/MutSel with codon bias         ###
### bias between Rate4Site vs dN/dS/MutSel with codon bias         ###
######################################################################

################### Correlations ###################
colfunc <- colorRampPalette(c("cyan2","navyblue"))
dNdS_p_bias_bl_cor_true <- ggplot(r_bias_dNdS,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
  ggtitle("dN/dS") +
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
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_bias_bl_cor_true <- ggplot(r_bias_MutSel,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
  ggtitle("MutSel") +
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
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))  

################### RMSD ###################
dNdS_p_bias_bl_rmsd_true <- ggplot(r_bias_dNdS,aes(bl,rmsd_true,colour=factor(num_taxa),group=num_taxa)) + 
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
  #scale_y_continuous(breaks=seq(0,1.5,0.25), limits = c(0,1.5))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_bias_bl_rmsd_true <- ggplot(r_bias_MutSel,aes(bl,rmsd_true,colour=factor(num_taxa),group=num_taxa)) + 
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
  #scale_y_continuous(breaks=seq(0,1.5,0.25), limits = c(0,1.5))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### bias ###################
dNdS_p_bias_bl_bias_true <- ggplot(r_bias_dNdS,aes(bl,bias_true,colour=factor(num_taxa),group=num_taxa)) + 
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
  ylab("Bias") +
  #scale_y_continuous(breaks=seq(-0.0010,0.0010,0.0002),
   #                  labels=c("-0.0010","-0.0008","-0.0006","-0.0004","-0.0002","0","0.0002","0.0004","0.0006","0.0008","0.0010"),limits = c(-0.0010,0.0010))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_bias_bl_bias_true <- ggplot(r_bias_MutSel,aes(bl,bias_true,colour=factor(num_taxa),group=num_taxa)) + 
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
  ylab("Bias") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(MutSel_p_bias_bl_bias_true)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_bias <- plot_grid(dNdS_p_bias_bl_cor_true+theme(legend.position="none",axis.title.x = element_blank()),
                  MutSel_p_bias_bl_cor_true+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
                  dNdS_p_bias_bl_rmsd_true+theme(legend.position="none",axis.title.x = element_blank()),
                  MutSel_p_bias_bl_rmsd_true+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
                  dNdS_p_bias_bl_bias_true+theme(legend.position="none"),
                  MutSel_p_bias_bl_bias_true+theme(legend.position="none",axis.title.y = element_blank()),
                  align = 'vh',
                  hjust = -1,
                  ncol=2,
                  nrow=3)

p_bias <- plot_grid( prow_bias, legend, rel_widths = c(2, .3))
p_bias

save_plot("plots/fig3_r4s_v_dNdS_MutSel_bias.png", p_bias,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

######################################################################
### Plotting Figure 4: Same as figure 2 but for inferred dN/dS     ###
### correlations between Rate4Site vs dN/dS/MutSel no codon bias   ###
### RMSD between Rate4Site vs dN/dS/MutSel no codon bias           ###
### bias between Rate4Site vs dN/dS/MutSel no codon bias           ###
######################################################################

################### Correlations ###################
colfunc <- colorRampPalette(c("cyan2","navyblue"))
dNdS_p_nobias_bl_cor_inferred <- ggplot(r_nobias_dNdS,aes(bl,cor_inferred,colour=factor(num_taxa),group=num_taxa)) + 
  ggtitle("dN/dS") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
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

MutSel_p_nobias_bl_cor_inferred <- ggplot(r_nobias_MutSel,aes(bl,cor_inferred,colour=factor(num_taxa),group=num_taxa)) + 
  ggtitle("MutSel") +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
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

################### RMSD ###################
dNdS_p_nobias_bl_rmsd_inferred <- ggplot(r_bias_dNdS,aes(bl,rmsd_inferred,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
               geom = "pointrange",
               size=0.6)+
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
  xlab("Branch Length") +
  ylab("RMSD") +
  scale_y_continuous(breaks=seq(0,1.5,0.25), limits = c(0,1.5))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_nobias_bl_rmsd_inferred <- ggplot(r_bias_MutSel,aes(bl,rmsd_inferred,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
               geom = "pointrange",
               size=0.6)+
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
  xlab("Branch Length") +
  ylab("RMSD") +
  #scale_y_continuous(breaks=seq(0,1.5,0.25), limits = c(0,1.5))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### bias ###################
dNdS_p_bias_bl_bias_true <- ggplot(r_bias_dNdS,aes(bl,bias,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
               geom = "pointrange",
               size=0.6)+
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
  xlab("Branch Length") +
  ylab("Bias") +
  scale_y_continuous(breaks=seq(-0.0010,0.0010,0.0002),
                     labels=c("-0.0010","-0.0008","-0.0006","-0.0004","-0.0002","0","0.0002","0.0004","0.0006","0.0008","0.0010"),limits = c(-0.0010,0.0010))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_bias_bl_nobias_inferred <- ggplot(r_nobias_MutSel,aes(bl,bias_inferred,colour=factor(num_taxa),group=num_taxa)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x)/sqrt(length(x)), 
               fun.ymax = function(x) mean(x) + sd(x)/sqrt(length(x)), 
               geom = "pointrange",
               size=0.6)+
  scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  scale_colour_manual(values=colfunc(5)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(num_taxa)),size=0.6)+
  xlab("Branch Length") +
  ylab("Bias") +
  #scale_y_continuous(breaks=seq(-0.0010,0.0010,0.0002),
   #                  labels=c("-0.0010","-0.0008","-0.0006","-0.0004","-0.0002","0","0.0002","0.0004","0.0006","0.0008","0.0010"),limits = c(-0.0010,0.0010))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(MutSel_p_bias_bl_bias_true)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_bias <- plot_grid(dNdS_p_bias_bl_cor_true+theme(legend.position="none",axis.title.x = element_blank()),
                       MutSel_p_bias_bl_cor_true+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
                       dNdS_p_bias_bl_rmsd_true+theme(legend.position="none",axis.title.x = element_blank()),
                       MutSel_p_bias_bl_rmsd_true+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_blank()),
                       dNdS_p_bias_bl_bias_true+theme(legend.position="none"),
                       MutSel_p_bias_bl_bias_true+theme(legend.position="none",axis.title.y = element_blank()),
                       align = 'vh',
                       hjust = -1,
                       ncol=2,
                       nrow=3)

p_bias <- plot_grid( prow_bias, legend, rel_widths = c(2, .3))
p_bias

save_plot("plots/fig2_r4s_v_dNdS_MutSel_bias.png", p_bias,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

#####################################################################
### Plotting                                                      ###
### correlations between Rate4Site vs dN/dS/MutSel no codon bias  ### 
### for number of taxa                                            ### 
#####################################################################

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

