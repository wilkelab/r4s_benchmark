library(dplyr)
library(tidyr)

t <- read.table("../test_sim_site_rates/site_rates.txt",header=T)
s <- read.table("../test_sim_site_rates/site_rates_info.txt",header=T)

ind <- match(t$Rate_Category,s$Rate_Category)

t$Rate_Probability <- s$Rate_Probability[ind]
t$Rate_Factor <- s$Rate_Factor[ind]
t <- t %>% separate(Rate_Factor, into = c("dN", "dS"), sep = "\\,") 
t <- select(t,-c(Partition_Index, Rate_Category))

write.table(t,"../test_sim_site_rates/site_rate_final.txt",sep="\t",quote=F,row.names=F)