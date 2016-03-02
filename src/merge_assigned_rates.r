library(dplyr)
library(tidyr)

model="mech_codon"
t1 <- list.files(paste0(model,"/sim_rates/assigned_rates/raw_rates"),pattern="sim_rates_rep",full.names=T)
t2 <- list.files(paste0(model,"/sim_rates/assigned_rates/raw_rates"),pattern="sim_rates_info_rep",full.names=T)
for (i in 1:length(t1)) {
  f1 <- t1[i] #sim_rates
  f2 <- t2[i] #sim_rates_info
  
  name_check1 <- substr(f1,60,120)
  name_check2 <- substr(f2,65,120)
  outfile <- paste0(model,"/sim_rates/assigned_rates/processed_rates/sim_rates_combined_rep",name_check1)
  if (name_check1==name_check2) {
    t <- read.table(f1,header=T)
    s <- read.table(f2,header=T)
    
    j <- inner_join(t,s,by="Rate_Category") %>%
     separate(Rate_Factor, into = c("dN", "dS"), sep = "\\,") %>%
      select(-c(Partition_Index.x,Rate_Category,Partition_Index.y,Model_Name))
    write.table(j,outfile,sep="\t",quote=F,row.names=F) 
  } else print("input file mismatch")
}