t <- read.table("../test_sim_site_rates/site_rates.txt",header=T)
s <- read.table("../test_sim_site_rates/site_rates_info.txt",header=T)

ind <- match(t$Rate_Category,s$Rate_Category)
rf_str <- as.character(s$Rate_Factor[ind])
dN <- 
dS <- 
#t$Rate_Probability <- s$Rate_Probability[ind]

write.table(t,"../test_sim_site_rates/site_rate_factor.txt",sep="\t",quote=F)