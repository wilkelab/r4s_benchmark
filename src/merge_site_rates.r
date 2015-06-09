library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
f1 <- as.character(args[1]) #site_rates
f2 <- as.character(args[2]) #site_rates_info
t <- read.table(f1,header=T)
s <- read.table(f2,header=T)

j <- inner_join(t,s,by="Rate_Category") %>%
  separate(Rate_Factor, into = c("dN", "dS"), sep = "\\,") %>%
  select(-c(Partition_Index.x,Rate_Category,Partition_Index.y,Model_Name))

write.table(t,f1,sep="\t",quote=F,row.names=F) #overwrite site_rate