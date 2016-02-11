library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

model="mech_codon"
r_bias <- read.csv(paste0(model,"/r4s_rates/processed_rates/all_r4s_norm_rates_bias.csv"),row.names=1)
r_nobias <- read.csv(paste0(model,"/r4s_rates/processed_rates/all_r4s_norm_rates_nobias.csv"),row.names=1)

# p1 <- ggplot(a,aes(dN,r4s_score)) + 
#   geom_point(size=1,alpha=0.7) + 
#   geom_smooth(method=lm) +
#   xlab(expression(bold("predicted rate (dN/dS)"))) +
#   ylab("rate4site score") +
#   theme(axis.text=element_text(size=8),legend.position="none") +
#   geom_text(aes(x=0.3,y=5,label=paste0(round(cor,2),sig),size=4)) +
#   scale_x_continuous(breaks=seq(0,1,0.5), labels=c("0","0.5","1.0"), limits = c(0.0,1.0)) + 
#   scale_y_continuous(breaks=seq(-2,6,2), limits = c(-2,6)) +
#   facet_grid(num_taxa ~ branch_len) +
#   background_grid(major = 'xy', minor = "none") + 
#   panel_border()
# ggsave(paste0("plots/",model,"_r4s_rates_v_sim_rates.png"))

##Plot Bias data
colfunc <- colorRampPalette(c("orange","darkred"))
p2 <- ggplot(r_bias,aes(num_taxa,cor,color=factor(bl),group=bl)) + 
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
ggsave(paste0("plots/",model,"_bias_cor_v_num_taxa.png"))

colfunc <- colorRampPalette(c("cyan2","navyblue"))
p3 <- ggplot(r_bias,aes(bl,cor,colour=factor(num_taxa),group=num_taxa)) + 
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
ggsave(paste0("plots/",model,"_bias_cor_v_branch_len.png"))

pg <- plot_grid(p2, p3, labels = c("A", "B"),ncol=2,nrow=1)
save_plot(paste0("plots/",model,"_bias_combined.png"), pg,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)

colfunc <- colorRampPalette(c("orange","darkred"))
p2 <- ggplot(r_nobias,aes(num_taxa,cor,color=factor(bl),group=bl)) + 
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
ggsave(paste0("plots/",model,"_nobias_cor_v_num_taxa.png"))

colfunc <- colorRampPalette(c("cyan2","navyblue"))
p3 <- ggplot(r_nobias,aes(bl,cor,colour=factor(num_taxa),group=num_taxa)) + 
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
ggsave(paste0("plots/",model,"_nobias_cor_v_branch_len.png"))

pg <- plot_grid(p2, p3, labels = c("A", "B"),ncol=2,nrow=1)
save_plot(paste0("plots/",model,"_nobias_combined.png"), pg,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 1, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
)


