library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

r_bias <- read.csv("other_methods_results/r4s_site_rates/processed_rates/all_r4s_orig_rates_bias.csv",row.names=1)
r_nobias <- read.csv("other_methods_results/r4s_site_rates/processed_rates/all_r4s_orig_rates_nobias.csv",row.names=1)

p1 <- ggplot(r_bias,aes(score)) + 
  geom_density()+
  scale_y_continuous(limit=c(0,4))+
  facet_grid(num_taxa~bl)+
  xlab("Rate4Site score") 

p2 <- ggplot(r_nobias,aes(score)) + 
  geom_density()+
  scale_y_continuous(limit=c(0,4))+
  facet_grid(num_taxa~bl)+
  xlab("Rate4Site score") 


r_bias$diff <- r_bias$score-r_bias$dnds
r_bias$scalar <- r_bias$dnds/r_bias$true2

p3 <- ggplot(r_bias,aes(scalar))+
  geom_histogram(binwidth=1)+
  facet_grid(ntaxa~bl)

