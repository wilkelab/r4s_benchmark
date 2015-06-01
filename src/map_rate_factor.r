library(dplyr)

t <- read.table("../test_sim_site_rates/site_rates.txt",header=T)
s <- read.table("../test_sim_site_rates/site_rates_info.txt",header=T)

ind <- match(t$Rate_Category,s$Rate_Category)
rf_str <- as.character(s$Rate_Factor[ind])
rf_v <- unlist(strsplit(rf_str,","))
dN <- rf_v[c(TRUE,FALSE)] 
dS <- rf_v[c(FALSE,TRUE)]

t$Rate_Probability <- s$Rate_Probability[ind]
t$dN <- dN
t$dS <- dS
t <- select(t,-c(Partition_Index, Rate_Category))

write.table(t,"../test_sim_site_rates/site_rate_final.txt",sep="\t",quote=F,row.names=F)