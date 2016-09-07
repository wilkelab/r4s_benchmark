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

prot_labels <- c(
  "hiv1_capsid" = "Capsid",
  "hiv1_gp120" = "gp120",
  "hiv1_integrase" = "Integrase",
  "hiv1_protease" = "Protease"
)
 
new_inf_lst <- rep(NA, length(r$inferred))
new_inf_lst[which(r$inferred < log(1.001))] = r$inferred[r$inferred < log(1.001)]
r$rug_inf <- new_inf_lst

new_score_lst <- rep(NA, length(r$inferred))
new_score_lst[which(r$inferred < log(1.001))] = r$score[r$inferred < log(1.001)]
r$rug_score <- new_score_lst

p <- ggplot(r,aes(score,inferred))+ 
  geom_point()+
  geom_rug(mapping=aes(rug_score,rug_inf),sides="b",na.rm=TRUE)+
 # geom_smooth(aes(score,inferred),method="lm",se=FALSE)+
  xlab("Rate4Site rates") +
  ylab("Inferred rates") +
  facet_wrap(~prot_name,labeller=as_labeller(prot_labels)) +
  guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
  geom_text(aes(x=0.01, y=100, label=paste0("rho=",round(cor,2))), size=4)+
  scale_x_log10(labels=fancy_scientific)+
  scale_y_log10(labels=fancy_scientific)+
  coord_cartesian(ylim=c(0.001, 100),xlim=c(0.003,100))+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))
print(p)
ggsave("plots/fig6_r4s_v_inf_rates_natural_prot.png",plot=p)
