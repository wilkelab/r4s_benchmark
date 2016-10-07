library(dplyr)
library(tidyr)

model="mech_codon"
t1 <- list.files(paste0(model,"/assigned_rates/raw_rates"),pattern="sim_rates_rep",full.names=T)

for (i in 1:length(t1)) {
  f1 <- t1[i] #sim_rates
  
  str <- regexpr("_\\w+.txt",f1)[1]
  end <- regexpr(".txt",f1)[1]
  bias <- substr(f1,str+1,end-1)
  
  str <- regexpr("bl\\d+",f1)[1]
  end <- regexpr("_\\w+.txt",f1)[1]
  bl <- as.numeric(substr(f1,str+2,end-1))

  str <- regexpr("rep\\d+",f1)[1]
  end <- regexpr("_n\\d+",f1)[1]
  rep <- as.numeric(substr(f1,str+3,end-1))

  str <- regexpr("n\\d+",f1)[1]
  end <- regexpr("_bl\\d+",f1)[1]
  n <- as.numeric(substr(f1,str+1,end-1))
  
  info_rate_file <- paste0("sim_rates_info_rep",rep,"_n",n,"_bl",bl,"_",bias,".txt")
  info_rate_t <- read.table(paste0(model,"/assigned_rates/raw_rates/",info_rate_file),header=T)
    
  rate_t <- read.table(f1,header=T)
  
  j <- rate_t %>% select(-c(Partition_Index)) %>% 
    inner_join(info_rate_t,by="Rate_Category") %>%
    separate(Rate_Factor, into = c("dN", "dS"), sep = "\\,") 

  outfile <- paste0(model,"/assigned_rates/processed_rates/sim_rates_combined_rep",rep,"_n",n,"_bl",bl,"_",bias,".txt")
  write.table(j,outfile,sep="\t",quote=F,row.names=F)
}