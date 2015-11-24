library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t_bias <- read.csv("other_methods_results/all_methods_processed_rates/all_methods_bias.csv")
t_nobias <- read.csv("other_methods_results/all_methods_processed_rates/all_methods_nobias.csv")

rf_bias <- t_bias %>% filter(method=="Rate4Site" | method == "FEL1")
rf_nobias <- t_nobias %>% filter(method=="Rate4Site" | method == "FEL1")

## plot number of taxa vs correlation
ggplot(rf_nobias,aes(ntaxa,cor1,color=factor(method),group=method)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange")+
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(method))) + 
  xlab("Number of Taxa") +
  ylab("Correation (Spearman)") +
  scale_x_log10(breaks=c(128,256,512,1024,2048)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) 

ggplot(rf_bias,aes(ntaxa,cor1,color=factor(method),group=method)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange")+
  stat_summary(fun.y = mean,geom = "line",aes(color=factor(method))) + 
  xlab("Number of Taxa") +
  ylab("Correation (Spearman)") +
  scale_x_log10(breaks=c(128,256,512,1024,2048)) +
  scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))

