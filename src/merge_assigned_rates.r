library(dplyr)
library(tidyr)

setwd("r4s_benchmark")
model="mech_codon"
distr_dir = 'gamma_distr'
#distr_dir = 'uniform_distr'

t1 <- list.files(paste0(model,"/assigned_rates/raw_rates/",distr_dir),pattern="sim_rates_rep",full.names=T)

for (i in 1:length(t1)) {
  f1 <- t1[i] #sim_rates

  info_rate_file <- gsub("sim_rates", "sim_rates_info", f1)
  info_rate_t <- read.table(info_rate_file,header=T)
  
  rate_t <- read.table(f1,header=T)
  
  j <- rate_t %>% select(-c(Partition_Index)) %>% 
    inner_join(info_rate_t,by="Rate_Category") %>%
    separate(Rate_Factor, into = c("dN", "dS"), sep = "\\,") 

  if (distr_dir=='gamma_distr') {
    outtemp <- gsub("sim_rates", "sim_gamma_rates_combined", f1)
    outfinal <- gsub("raw_rates/gamma_distr", "processed_rates", outtemp)
    write.table(j,outfinal,sep="\t",quote=F,row.names=F)
  } else {
    outtemp <- gsub("sim_rates", "sim_rates_combined", f1)
    outfinal <- gsub("raw_rates/uniform_distr", "processed_rates", outtemp)
    write.table(j,outfinal,sep="\t",quote=F,row.names=F)
  }
}