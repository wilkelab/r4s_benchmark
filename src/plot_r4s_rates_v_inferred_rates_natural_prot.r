library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(readr)

setwd("r4s_benchmark/")
r <- read_csv("natural_prot/processed_rates/all_orig_rates.csv")

r1 <- r %>% group_by(prot_name) %>% mutate(inferred_norm=inferred/mean(inferred)) 

##extract inferred rates and scores to be plotted in the geom_rug() 
new_inf_lst <- rep(NA, length(r$inferred))
new_inf_lst[which(r$inferred < log(1.001))] = r$inferred[r$inferred < log(1.001)]
r$rug_inf <- new_inf_lst

new_score_lst <- rep(NA, length(r$inferred))
new_score_lst[which(r$inferred < log(1.001))] = r$score[r$inferred < log(1.001)]
r$rug_score <- new_score_lst

###############################################################################
### Plot sub figure 6:                                                      ###
### correlations between Rate4Site score vs inferred dN/dS for HIV 1 Capsid ###
###############################################################################
r_capsid <- filter(r1,prot_name=="hiv1_capsid")

c <- cor.test(r_capsid$score, r_capsid$inferred, method='sp')
label <- paste0("ρ = ", signif(c$estimate, 2))

p_capsid <- ggplot(r_capsid,aes(score,inferred_norm))+ 
  geom_point()+
  ggtitle("Capsid")+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
 # geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  scale_x_log10(breaks=c(0.01,0.1,1,10),label=c("0.01","0.1","1","10"))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.001, 13),xlim=c(0.003,10))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title.x=element_blank())
p1 <- p_capsid + draw_text(label, 0.05, 0.3)

###############################################################################
### Plot sub figure 6:                                                      ###
### correlations between Rate4Site score vs inferred dN/dS for HIV 1 gp120  ###
###############################################################################
r_gp120 <- filter(r1,prot_name=="hiv1_gp120")

c <- cor.test(r_gp120$score, r_gp120$inferred, method='sp')
label <- paste0("ρ = ", signif(c$estimate, 2))

p_gp120 <- ggplot(r_gp120,aes(score,inferred_norm))+ 
  geom_point()+
  ggtitle("gp120")+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
  # geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  scale_x_log10(breaks=c(0.01,0.1,1,10),label=c("0.01","0.1","1","10"))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.001, 13),xlim=c(0.003,10))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
p2 <- p_gp120 + draw_text(label, 0.05, 0.3)

##################################################################################
### Plot sub figure 6:                                                         ###
### correlations between Rate4Site score vs inferred dN/dS for HIV 1 integrase ###
##################################################################################
r_integrase <- filter(r1,prot_name=="hiv1_integrase")

c <- cor.test(r_integrase$score, r_integrase$inferred, method='sp')
label <- paste0("ρ = ", signif(c$estimate, 2))

p_gp120 <- ggplot(r_gp120,aes(score,inferred_norm))+ 
  geom_point()+
  ggtitle("Integrase")+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
  # geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  scale_x_log10(breaks=c(0.01,0.1,1,10),label=c("0.01","0.1","1","10"))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.001, 13),xlim=c(0.003,10))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
p3 <- p_gp120 + draw_text(label, 0.05, 0.3)

#################################################################################
### Plot sub figure 6:                                                        ###
### correlations between Rate4Site score vs inferred dN/dS for HIV 1 protease ###
#################################################################################
r_protease <- filter(r1,prot_name=="hiv1_protease")

c <- cor.test(r_protease$score, r_protease$inferred, method='sp')
label <- paste0("ρ = ", signif(c$estimate, 2))

p_protease <- ggplot(r_protease,aes(score,inferred_norm))+ 
  geom_point()+
  ggtitle("Protease")+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  #geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
  # geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  scale_x_log10(breaks=c(0.01,0.1,1,10),label=c("0.01","0.1","1","10"))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.001, 13),xlim=c(0.003,10))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title.y=element_blank())
p4 <- p_protease + draw_text(label, 0.05, 0.3)

p_final <- plot_grid(p1,
                         p2,
                         p3,
                         p4,
                         labels = "AUTO",
                         align = 'vh',
                         hjust = -1,
                         ncol=2,
                         nrow=2)
save_plot("plots/fig6_r4s_v_inf_rates_natural_prot.png", p_final,
          ncol = 2, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3)

