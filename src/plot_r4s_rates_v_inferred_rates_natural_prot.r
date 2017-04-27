library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(readr)

setwd("r4s_benchmark/")
r <- read_csv("natural_prot/processed_rates/all_orig_rates.csv")

name_df <- data.frame(prot_name=c("ENST00000000412",
                                  "ENST00000009530",
                                  "ENST00000011653",
                                  "ENST00000014914",
                                  "ENST00000023897",
                                  "ENST00000053243",
                                  "hiv1_capsid",
                                  "hiv1_gp120",
                                  "hiv1_integrase",
                                  "hiv1_matrix",
                                  "hiv1_protease",
                                  "hiv1_rt"),
                      aln_length=c(279, #ENST00000000412
                                   305, #ENST00000009530
                                   546, #ENST00000011653
                                   367, #ENST00000014914
                                   469, #ENST00000023897
                                   199, #ENST00000053243
                                   250, #hiv1_capsid
                                   750, #hiv1_gp120
                                   328, #hiv1_integrase
                                   172, #hiv1_matrix
                                   106,#hiv1_protease
                                   629), #hiv1_rt
                      plot_prot_name=c("M6PR", ##Mannose-6-phosphate receptor
                                       "CD74",
                                       "CD4",
                                       "GPRC5A", #full name "G protein-coupled\nreceptor class C"
                                       "GABRA1", ##Gamma-aminobutyric acid type A receptor
                                       "TNFRSF17",
                                       "Capsid",
                                       "gp120",
                                       "Integrase",
                                       "Matrix",
                                       "Protease",
                                       "Reverse Transcriptase"))

r1 <- r %>% 
  inner_join(name_df) %>%
  group_by(prot_name) %>%
  mutate(score_norm=score/mean(score),
         panel_label=paste0(plot_prot_name[1],"\nn = ",num_taxa[1],", l = ",n(pos)))

r1 <- r1 %>% 
  na.omit()  %>%
  summarize(cor=cor(inferred,score_norm,method="spearman",use="pairwise.complete.obs"),
         p_value=cor.test(inferred,score)$p.value,
         cor_label=paste("ρ =", round(cor, 2))) %>% 
  select(cor_label,prot_name) %>%
  left_join(r1)


##extract inferred rates and scores to be plotted in the geom_rug() 
new_inf_lst <- rep(NA, length(r1$inferred))
new_inf_lst[which(r1$inferred < .001)] = r1$inferred[r1$inferred < 0.001]
r1$rug_inf <- new_inf_lst
 
filter_inf <- r1$inferred
filter_inf[which(r1$inferred < .001)] = rep(NA, length(which(r1$inferred < .001)))
r1$inferred_filtered <- filter_inf
   
new_score_lst <- rep(NA, length(r1$inferred))
new_score_lst[which(r1$inferred < 0.001)] = r1$score_norm[r1$inferred < 0.001]
r1$rug_score <- new_score_lst

filter_score_norm <- r1$score_norm
filter_score_norm [which(r1$inferred < .001)] = rep(NA, length(which(r1$inferred < .001)))
r1$score_filtered_norm <- filter_score_norm 

##################################################################################
### Plot figure 6:                                                             ###
### correlations between Rate4Site score vs inferred dN/dS for HIV 1 integrase ###
##################################################################################
r1_gpcr <- r1 %>% filter(prot_name=="ENST00000000412" |
                         prot_name=="ENST00000009530" | 
                         prot_name=="ENST00000011653" |
                         prot_name=="ENST00000014914" |
                         prot_name=="ENST00000023897" |
                          prot_name=="ENST00000053243")
 
 # r1_gpcr_cor <- r_cor %>% filter(prot_name=="ENST00000000412" |
 #                            prot_name=="ENST00000009530" | 
 #                            prot_name=="ENST00000011653" |
 #                            prot_name=="ENST00000014914" |
 #                            prot_name=="ENST00000023897" |
 #                              prot_name=="ENST00000053243")
  
r1_hiv <- r1 %>% filter(prot_name=="hiv1_capsid" |
                        prot_name=="hiv1_gp120" |
                        prot_name=="hiv1_integrase" |
                        prot_name=="hiv1_matrix" |
                        prot_name=="hiv1_protease" |
                        prot_name=="hiv1_rt")

# r1_hiv_cor <- r_cor %>% filter(prot_name=="hiv1_capsid" |
#                           prot_name=="hiv1_gp120" |
#                           prot_name=="hiv1_integrase" |
#                           prot_name=="hiv1_matrix" |
#                           prot_name=="hiv1_protease" |
#                           prot_name=="hiv1_rt")

p_gpcr <- ggplot(r1_gpcr,aes(score_filtered_norm,inferred_filtered))+ 
  geom_point()+
  #geom_hline(yintercept=1)+
  #ggtitle(prot_name)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  facet_wrap(~panel_label)+
  scale_x_log10(breaks=c(0.1,1,10),label=c("0.1","1","10"))+
  scale_y_log10(breaks=c(0.01,0.1,1,10),label=c("0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.01, 13),xlim=c(0.09,10))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))+
  geom_text(aes(0.21,9,label=cor_label))

p_hiv <- ggplot(r1_hiv,aes(score_filtered_norm,inferred_filtered))+ 
  geom_point()+
  #ggtitle(prot_name)+
  geom_abline(slope=1,intercept=0)+
  #background_grid("xy")+
  geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
  #geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site scores") +
  ylab("Inferred dN/dS") +
  facet_wrap(~panel_label)+
  scale_x_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  scale_y_log10(breaks=c(0.001,0.01,0.1,1,10),label=c("0.001","0.01","0.1","1","10"))+
  coord_cartesian(ylim=c(0.0007, 13),xlim=c(0.0007,16))+
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 10))+
  geom_text(aes(0.004,9,label=cor_label))

prow <- plot_grid(p_gpcr+ggtitle("Membrane proteins"),
                       p_hiv+ggtitle("HIV 1 proteins"),
                       labels=c("A","B"),
                       align = 'vh',
                       hjust = -1,
                       ncol=1,
                       nrow=2)

save_plot("plots/r4s_v_dNdS_inferred_natural_prot.pdf", prow,
          ncol = 1, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_height=5,
         base_width=7)
