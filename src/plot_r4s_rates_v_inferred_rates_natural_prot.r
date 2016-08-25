library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(readr)

setwd("r4s_benchmark/")
r <- read_csv("natural_prot/processed_rates/all_orig_rates.csv")

#####################################################################
### Calculating correlations, RMSD, and bias                      ###
#####################################################################

#norm_r <- r %>% group_by(prot_name) %>% mutate(mean_score=mean(score),mean_inferred=mean(inferred)) %>% mutate(score_norm=score/mean_score,inferred_norm=inferred/mean_inferred)
r <- r %>% group_by(prot_name) %>% mutate(cor=cor(score,inferred,method="spearman",use="pairwise.complete.obs"))

p <- ggplot(r,aes(score,inferred))+ 
  geom_point()+
  #scale_x_log10(breaks=c(0.0025,0.01,0.04,0.16,0.64),labels=c("0.0025","0.01","0.04","0.16","0.64")) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  xlab("Rate4Site rates") +
  ylab("Inferred rates") +
  facet_wrap(~prot_name) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_text(aes(x=1, y=15, label=paste0("rho=",round(cor,2))), size=6)+
  #coord_cartesian(ylim=c(0, 1),xlim=c(0.0023,0.66))+
  #scale_y_continuous(breaks=seq(0,1,0.2))+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

ggsave("plots/fig6_r4s_v_inf_rates_natural_prot.png",plot=p)
