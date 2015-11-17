library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

model_lst = c("mech_codon_dN","mech_codon_dN_dS","mut_sel_dN","mut_sel_dN_dS")

##for (model in model_lst) {
t1 <- list.files(paste0(model,"/r4s_site_rates/norm_rates"),full.names=T)
info = file.info(t1)
t1 <- t1[info$size != 0]

if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
  t2 <- list.files(paste0(model,"/sim_site_rates/merged_output"),
                   full.names=T)
}
