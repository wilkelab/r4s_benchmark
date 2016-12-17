library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(readr)

setwd("r4s_benchmark/")
r <- read_csv("natural_prot/processed_rates/all_orig_rates.csv")

name_df <- data.frame(prot_name=unique(r$prot_name),
                      plot_prot_name=c("Mannose-6-ph. rec.", ##Mannose-6-phosphate receptor
                                       "CD74",
                                       "CD4",
                                       "G prot.-coupled rec.", #full name "G protein-coupled\nreceptor class C"
                                       "Gamma-aminobut. acid rec.", ##Gamma-aminobutyric acid type A receptor
                                       "Capsid",
                                       "gp120",
                                       "Integrase",
                                       "Matrix",
                                       "Protease",
                                       "Reverse Transcriptase"))

r1 <- r %>% inner_join(name_df)
  
r_cor <- r1 %>% 
  na.omit() %>%
  group_by(plot_prot_name,prot_name) %>%
  summarize(cor=cor(inferred,score,method="spearman",use="pairwise.complete.obs"),
         p_value=cor.test(inferred,score)$p.value,
         cor_label=paste("ρ =", round(cor, 2))) 

##extract inferred rates and scores to be plotted in the geom_rug() 
 new_inf_lst <- rep(NA, length(r1$inferred))
 new_inf_lst[which(r1$inferred < .001)] = r1$inferred[r1$inferred < 0.001]
 r1$rug_inf <- new_inf_lst
 
 filter_inf <- r1$inferred
 filter_inf[which(r1$inferred < .001)] = rep(NA, length(which(r1$inferred < .001)))
 r1$inferred_filtered <- filter_inf
   
new_score_lst <- rep(NA, length(r1$inferred))
new_score_lst[which(r1$inferred < 0.001)] = r1$score[r1$inferred < 0.001]
r1$rug_score <- new_score_lst

filter_score <- r1$score
filter_score[which(r1$inferred < .001)] = rep(NA, length(which(r1$inferred < .001)))
r1$score_filtered <- filter_score

##################################################################################
### Plot figure 6:                                                             ###
### correlations between Rate4Site score vs inferred dN/dS for HIV 1 integrase ###
##################################################################################
r1_gpcr <- r1 %>% filter(prot_name=="ENST00000000412" |
                         prot_name=="ENST00000009530" | 
                         prot_name=="ENST00000011653" |
                         prot_name=="ENST00000014914" |
                         prot_name=="ENST00000023897")
 
 r1_gpcr_cor <- r_cor %>% filter(prot_name=="ENST00000000412" |
                            prot_name=="ENST00000009530" | 
                            prot_name=="ENST00000011653" |
                            prot_name=="ENST00000014914" |
                            prot_name=="ENST00000023897")

r1_hiv <- r1 %>% filter(prot_name=="hiv1_capsid" |
                        prot_name=="hiv1_gp120" |
                        prot_name=="hiv1_integrase" |
                        prot_name=="hiv1_matrix" |
                        prot_name=="hiv1_protease" |
                        prot_name=="hiv1_rt")

r1_hiv_cor <- r_cor %>% filter(prot_name=="hiv1_capsid" |
                          prot_name=="hiv1_gp120" |
                          prot_name=="hiv1_integrase" |
                          prot_name=="hiv1_matrix" |
                          prot_name=="hiv1_protease" |
                          prot_name=="hiv1_rt")

p_gpcr <- ggplot(r1_gpcr,aes(score_filtered,inferred_filtered))+ 
  geom_point()+
  #geom_hline(yintercept=1)+
  #ggtitle(prot_name)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  facet_wrap(~plot_prot_name)+
  scale_x_log10(breaks=c(0.1,1,10),label=c("0.1","1","10"))+
  scale_y_log10(breaks=c(0.01,0.1,1,10),label=c("0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.01, 13),xlim=c(0.09,12))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))+
  geom_text(data=r1_gpcr_cor,aes(0.21,9,label=cor_label), inherit.aes=FALSE)

p_hiv <- ggplot(r1_hiv,aes(score_filtered,inferred_filtered))+ 
  geom_point()+
  #ggtitle(prot_name)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  facet_wrap(~plot_prot_name)+
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.0007, 13),xlim=c(0.0007,12))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))+
  geom_text(data=r1_hiv_cor,aes(0.004,9,label=cor_label), inherit.aes=FALSE)

prow <- plot_grid(p_gpcr+ggtitle("Membrane proteins"),
                       p_hiv+ggtitle("HIV 1 proteins"),
                       labels=c("A","B"),
                       align = 'vh',
                       hjust = -1,
                       ncol=1,
                       nrow=2)

save_plot("plots/fig6_r4s_v_inf_rates_natural_prot.png", prow,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_height=5,
         base_width=7)
