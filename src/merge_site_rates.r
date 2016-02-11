library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
f1 <- as.character(args[1]) #sim_rates
f2 <- as.character(args[2]) #sim_rates_info
out_f <- as.character(args[3])
t <- read.table(f1,header=T)
s <- read.table(f2,header=T)

j <- inner_join(t,s,by="Rate_Category") %>%
  separate(Rate_Factor, into = c("dN", "dS"), sep = "\\,") %>%
  select(-c(Partition_Index.x,Rate_Category,Partition_Index.y,Model_Name))

write.table(j,out_f,sep="\t",quote=F,row.names=F) #overwrite site_rate