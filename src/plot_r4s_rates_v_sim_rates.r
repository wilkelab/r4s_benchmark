library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

r_bias_dNdS <- read.csv("mech_codon/r4s_rates/processed_rates/all_r4s_orig_rates_bias.csv",row.names=1)
r_nobias_dNdS <- read.csv("mech_codon/r4s_rates/processed_rates/all_r4s_orig_rates_nobias.csv",row.names=1)
r_bias_MutSel <- read.csv("mut_sel/r4s_rates/processed_rates/all_r4s_orig_rates_bias.csv",row.names=1)
r_nobias_MutSel <- read.csv("mut_sel/r4s_rates/processed_rates/all_r4s_orig_rates_nobias.csv",row.names=1)

#####################################################################
### Plotting Figure 2:                                            ###
### correlations between Rate4Site vs dN/dS/MutSel no codon bias  ###
### RMSD between Rate4Site vs dN/dS/MutSel no codon bias          ###
### bias between Rate4Site vs dN/dS/MutSel no codon bias          ###
#####################################################################

################### Correlations ###################
colfunc <- colorRampPalette(c("cyan2","navyblue"))
dNdS_p_nobias_bl_cor_true <- ggplot(r_nobias_dNdS,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
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

MutSel_p_nobias_bl_cor_true <- ggplot(r_nobias_MutSel,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
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
dNdS_p_nobias_bl_rmsd_true <- ggplot(r_nobias_dNdS,aes(bl,rmsd,colour=factor(num_taxa),group=num_taxa)) + 
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
  scale_y_continuous(breaks=seq(0,0.6,0.1), limits = c(0,0.6))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_nobias_bl_rmsd_true <- ggplot(r_nobias_MutSel,aes(bl,rmsd,colour=factor(num_taxa),group=num_taxa)) + 
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
  scale_y_continuous(breaks=seq(0,0.6,0.1), limits = c(0,0.6))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### bias ###################
dNdS_p_nobias_bl_bias_true <- ggplot(r_nobias_dNdS,aes(bl,bias,colour=factor(num_taxa),group=num_taxa)) + 
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
  scale_y_continuous(breaks=seq(-0.0008,0.0008,0.0002),labels=c("-0.0008","-0.0006","-0.0004","-0.0002","0","0.0002","0.0004","0.0006","0.0008"),limits = c(-0.0008,0.0008))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_nobias_bl_bias_true <- ggplot(r_nobias_MutSel,aes(bl,bias,colour=factor(num_taxa),group=num_taxa)) + 
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
  scale_y_continuous(breaks=seq(-0.0008,0.0008,0.0002),labels=c("-0.0008","-0.0006","-0.0004","-0.0002","0","0.0002","0.0004","0.0006","0.0008"),limits = c(-0.0008,0.0008))+ 
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(MutSel_p_nobias_bl_bias_true)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow <- plot_grid(dNdS_p_nobias_bl_cor_true+theme(legend.position="none"),
                MutSel_p_nobias_bl_cor_true+theme(legend.position="none"),
                dNdS_p_nobias_bl_rmsd_true+theme(legend.position="none"),
                MutSel_p_nobias_bl_rmsd_true+theme(legend.position="none"),
                dNdS_p_nobias_bl_bias_true+theme(legend.position="none"),
                MutSel_p_nobias_bl_bias_true+theme(legend.position="none"),
                align = 'vh',
                hjust = -1,
                ncol=2,
                nrow=3)
prow
p <- plot_grid( prow, legend, rel_widths = c(3, .3))

save_plot("plots/fig2_r4s_v_dNdS_MutSel_nobias.png", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 3, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)



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

