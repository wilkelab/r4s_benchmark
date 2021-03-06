library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(readr)

setwd("r4s_benchmark/")
r_bias_dNdS_raw <- read_csv("mech_codon/processed_rates/all_r4s_orig_rates_bias.csv")
r_nobias_dNdS_raw <- read_csv("mech_codon/processed_rates/all_r4s_orig_rates_nobias.csv")
r_bias_MutSel_raw <- read_csv("mut_sel/processed_rates/all_r4s_orig_rates_bias.csv")
r_nobias_MutSel_raw <- read_csv("mut_sel/processed_rates/all_r4s_orig_rates_nobias.csv")

#####################################################################
### Calculating correlations, RMSD, and bias                      ###
#####################################################################


r_bias_dNdS <- r_bias_dNdS_raw %>% 
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>% 
  summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
            rmsd_true=sqrt(mean((score_norm - true_norm)^2)),
            cor_inferred=cor(score,inferred,method="spearman",use="pairwise.complete.obs"),
            rmsd_inferred=sqrt(mean((score_norm - inferred_norm)^2)))

r_nobias_dNdS <- r_nobias_dNdS_raw %>%
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>% 
  summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
            rmsd_true=sqrt(mean((score_norm - true_norm)^2)),
            cor_inferred=cor(score,inferred,method="spearman",use="pairwise.complete.obs"),
            rmsd_inferred=sqrt(mean((score_norm - inferred_norm)^2)))

r_bias_MutSel <- r_bias_MutSel_raw %>% 
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>% 
  summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
            rmsd_true=sqrt(mean((score_norm - true_norm)^2)),
            cor_inferred=cor(score,inferred,method="spearman",use="pairwise.complete.obs"),
            rmsd_inferred=sqrt(mean((score_norm - inferred_norm)^2)))

r_nobias_MutSel <- r_nobias_MutSel_raw %>% 
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>% 
  summarise(cor_true=cor(score,true,method="spearman",use="pairwise.complete.obs"),
            rmsd_true=sqrt(mean((score_norm - true_norm)^2)),
            cor_inferred=cor(score,inferred,method="spearman",use="pairwise.complete.obs"),
            rmsd_inferred=sqrt(mean((score_norm - inferred_norm)^2)))

rr_nobias_dNdS <- r_nobias_dNdS_raw %>% 
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>%
  filter(rep==1,num_taxa==512,bl==0.0025 | bl==0.04 | bl==0.64)

rr_bias_dNdS <- r_bias_dNdS_raw %>% 
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>%
  filter(rep==1,num_taxa==512,bl==0.0025 | bl==0.04 | bl==0.64)

rr_nobias_MutSel <- r_nobias_MutSel_raw %>% 
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>%
  filter(rep==1,num_taxa==512,bl==0.0025 | bl==0.04 | bl==0.64)

rr_bias_MutSel <- r_bias_MutSel_raw %>% 
  na.omit() %>%
  group_by(bl,num_taxa,rep) %>% 
  mutate(score_norm=score/mean(score),true_norm=true/mean(true),inferred_norm=inferred/mean(inferred)) %>%
  filter(rep==1,num_taxa==512,bl==0.0025 | bl==0.04 | bl==0.64)

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

#########################################################################
### Plotting:                                                         ###
### Site-wise Rate4Site scores vs site-wise dN/dS for the dN/dS model ###
#########################################################################

### Site-wise normalized R4S vs normalized dN/dS ####
r4s_vs_dNdS_nobias <- ggplot(rr_nobias_dNdS,aes(true_norm,score_norm))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("True dN/dS (normalized)") +
  ylab("Rate4Site scores") +
  facet_wrap(~bl)+
  scale_x_continuous(breaks=seq(0,2,0.5),
                     label=c("0","0.5","1.0","1.5","2.0"),
                     expand = c(0.02, 0.02))+
  scale_y_continuous(breaks=seq(0,3.5,0.5),
                     label=c("0","0.5","1.0","1.5","2.0","2.5","3.0","3.5"),
                     expand = c(0.02, 0.02))+
  coord_cartesian(xlim=c(0,2),ylim=c(0,3.51))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

r4s_vs_dNdS_bias <- ggplot(rr_bias_dNdS,aes(true_norm,score_norm))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("True dN/dS (normalized)") +
  ylab("Rate4Site scores") +
  facet_wrap(~bl)+
  scale_x_continuous(breaks=seq(0,2,0.5),
                     label=c("0","0.5","1.0","1.5","2.0"),
                     expand = c(0.02, 0.02))+
  scale_y_continuous(breaks=seq(0,3.5,0.5),
                     label=c("0","0.5","1.0","1.5","2.0","2.5","3.0","3.5"),
                     expand = c(0.02, 0.02))+
  coord_cartesian(xlim=c(0,2),ylim=c(0,3.51))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

p <- plot_grid(r4s_vs_dNdS_nobias,
               r4s_vs_dNdS_bias,
               labels="AUTO",
               align = 'vh',
               hjust = -1,
               ncol=1,
               nrow=2)

save_plot("plots/site-wise_r4s_v_dNdS_true.pdf", p,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_height=4,
          base_width=10)

#########################################################################
### Plotting:                                                         ###
### Site-wise Rate4Site scores vs site-wise dN/dS for the MutSel model###
#########################################################################

### Site-wise normalized R4S vs normalized dN/dS ####
r4s_vs_MutSel_nobias <- ggplot(rr_nobias_MutSel,aes(true_norm,score_norm))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("True dN/dS (normalized)") +
  ylab("Rate4Site scores") +
  facet_wrap(~bl)+
  scale_x_continuous(breaks=seq(0,2,0.5),
                     label=c("0","0.5","1.0","1.5","2.0"),
                     expand = c(0.02, 0.02))+
  scale_y_continuous(breaks=seq(0,4,1),
                     label=c("0","1.0","2.0","3.0","4.0"),
                     expand = c(0.02, 0.02))+
  coord_cartesian(xlim=c(0,2),ylim=c(0,4.1))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

r4s_vs_MutSel_bias <- ggplot(rr_bias_MutSel,aes(true_norm,score_norm))+
  geom_point(size=2)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("True dN/dS (normalized)") +
  ylab("Rate4Site scores") +
  facet_wrap(~bl)+
  scale_x_continuous(breaks=seq(0,2,0.5),
                     label=c("0","0.5","1.0","1.5","2.0"),
                     expand = c(0.02, 0.02))+
  scale_y_continuous(breaks=seq(0,4,1),
                     label=c("0","1.0","2.0","3.0","4.0"),
                     expand = c(0.02, 0.02))+
  coord_cartesian(xlim=c(0,2),ylim=c(0,4.1))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

p <- plot_grid(r4s_vs_MutSel_nobias,
               r4s_vs_MutSel_bias,
               labels="AUTO",
               align = 'vh',
               hjust = -1,
               ncol=1,
               nrow=2)

save_plot("plots/site-wise_r4s_v_MutSel_true.pdf", p,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_height=4,
          base_width=10)

#########################################################################
### Plotting:                                                         ###
### correlations between Rate4Site vs dN/dS no codon bias/codon bias  ###
### RMSD between Rate4Site vs dN/dS no codon bias/codon bias          ###
#########################################################################

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

dNdS_p_bias_bl_cor_true <- ggplot(r_bias_dNdS,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
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
  coord_cartesian(ylim=c(0,0.8),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,0.8,0.2))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))
        
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
  coord_cartesian(ylim=c(0,0.8),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,0.8,0.2))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(dNdS_p_nobias_bl_cor_true)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_bias <- plot_grid(dNdS_p_nobias_bl_cor_true+theme(legend.position="none")+ggtitle("Constant dS"),
                       dNdS_p_bias_bl_cor_true+theme(legend.position="none",axis.title.y = element_blank())+ggtitle("Variable dS"),
                       dNdS_p_nobias_bl_rmsd_true+theme(legend.position="none"),
                       dNdS_p_bias_bl_rmsd_true+theme(legend.position="none",axis.title.y = element_blank()),
                       labels=c("A","C","B","D"),
                       align = 'vh',
                       hjust = -1,
                       ncol=2,
                       nrow=2)

p <- plot_grid( prow_bias, legend, rel_widths = c(2, .3))

save_plot("plots/r4s_v_dNdS_true.pdf", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

#########################################################################
### Plotting Figure 3:                                                ###
### correlations between Rate4Site vs MutSel no codon bias/codon bias ###
### RMSD between Rate4Site vs MutSel no codon bias/codon bias         ###
#########################################################################

############ Correlations vs Branch len ###########
colfunc <- colorRampPalette(c("cyan2","navyblue"))
        
MutSel_p_nobias_bl_cor_true <- ggplot(r_nobias_MutSel,aes(bl,cor_true,colour=factor(num_taxa)))+ 
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


MutSel_p_bias_bl_cor_true <- ggplot(r_bias_MutSel,aes(bl,cor_true,colour=factor(num_taxa),group=num_taxa)) + 
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
MutSel_p_nobias_bl_rmsd_true <- ggplot(r_nobias_MutSel,aes(bl,rmsd_true,colour=factor(num_taxa)))+ 
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
  coord_cartesian(ylim=c(0,1),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
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
  coord_cartesian(ylim=c(0,1),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))


################### Plotting the grid ###################
grobs <- ggplotGrob(MutSel_p_bias_bl_rmsd_true)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_bias <- plot_grid(MutSel_p_nobias_bl_cor_true+theme(legend.position="none")+ggtitle("Neutral synonymous codons"),
                       MutSel_p_bias_bl_cor_true+theme(legend.position="none",axis.title.y = element_blank())+ggtitle("Non-neutral synonymous codons"),
                       MutSel_p_nobias_bl_rmsd_true+theme(legend.position="none"),
                       MutSel_p_bias_bl_rmsd_true+theme(legend.position="none",axis.title.y = element_blank()),
                       labels="AUTO",
                       align = 'vh',
                       hjust = -1,
                       ncol=2,
                       nrow=2)

p <- plot_grid(prow_bias, legend, rel_widths = c(2, .3))

save_plot("plots/r4s_v_MutSel_true.pdf", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

########################################################################
### Plotting Figure 4: Same as figure 2 but for inferred dN/dS       ###
### correlations between Rate4Site vs dN/dS no codon bias/codon bias ###
### RMSD between Rate4Site vs dN/dS no codon bias/codon bias         ###
### bias between Rate4Site vs dN/dS no codon bias/codon bias         ###
########################################################################

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

dNdS_p_bias_bl_cor_inferred <- ggplot(r_bias_dNdS,aes(bl,cor_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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
  coord_cartesian(ylim=c(0,1.5),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,1.5,0.25))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

dNdS_p_bias_bl_rmsd_inferred <- ggplot(r_bias_dNdS,aes(bl,rmsd_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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
  coord_cartesian(ylim=c(0,1.5),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,1.5,0.25))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(dNdS_p_bias_bl_rmsd_inferred)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_bias <- plot_grid(dNdS_p_nobias_bl_cor_inferred+theme(legend.position="none")+ggtitle("Constant dS"),
                       dNdS_p_bias_bl_cor_inferred+theme(legend.position="none",axis.title.y = element_blank())+ggtitle("Variable dS"),
                       dNdS_p_nobias_bl_rmsd_inferred+theme(legend.position="none"),
                       dNdS_p_bias_bl_rmsd_inferred+theme(legend.position="none",axis.title.y = element_blank()),
                       labels=c("A","C","B","D"),
                       align = 'vh',
                       hjust = -1,
                       ncol=2,
                       nrow=2)

p <- plot_grid( prow_bias, legend, rel_widths = c(2, .3))

save_plot("plots/r4s_v_dNdS_inferred.pdf", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

#########################################################################
### Plotting Figure 5: Same as figure 2 but for inferred dN/dS        ###
### correlations between Rate4Site vs MutSel no codon bias/codon bias ###
### RMSD between Rate4Site vs MutSel no codon bias/codon bias         ###
### bias between Rate4Site vs MutSel no codon bias/codon bias         ###
#########################################################################

############ Correlations vs Branch len ###########
colfunc <- colorRampPalette(c("cyan2","navyblue"))
MutSel_p_nobias_bl_cor_inferred <- ggplot(r_nobias_MutSel,aes(bl,cor_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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

MutSel_p_bias_bl_cor_inferred <- ggplot(r_bias_MutSel,aes(bl,cor_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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
MutSel_p_nobias_bl_rmsd_inferred <- ggplot(r_nobias_MutSel,aes(bl,rmsd_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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
  coord_cartesian(ylim=c(0,2),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,2,0.5))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

MutSel_p_bias_bl_rmsd_inferred <- ggplot(r_bias_MutSel,aes(bl,rmsd_inferred,colour=factor(num_taxa),group=num_taxa)) + 
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
  coord_cartesian(ylim=c(0,2),xlim=c(0.0023,0.66))+
  scale_y_continuous(breaks=seq(0,2,0.5))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

################### Plotting the grid ###################
grobs <- ggplotGrob(MutSel_p_bias_bl_rmsd_inferred )$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

prow_bias <- plot_grid(MutSel_p_nobias_bl_cor_inferred+theme(legend.position="none")+ggtitle("Neutral synonymous codons"),
                       MutSel_p_bias_bl_cor_inferred+theme(legend.position="none",axis.title.y = element_blank())+ggtitle("Non-neutral synonymous codons"),
                       MutSel_p_nobias_bl_rmsd_inferred+theme(legend.position="none"),
                       MutSel_p_bias_bl_rmsd_inferred+theme(legend.position="none",axis.title.y = element_blank()),
                       labels="AUTO",
                       align = 'vh',
                       hjust = -1,
                       ncol=2,
                       nrow=2)

p <- plot_grid( prow_bias, legend, rel_widths = c(2, .3))

save_plot("plots/r4s_v_MutSel_inferred.pdf", p,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)