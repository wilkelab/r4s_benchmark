library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
t <- read.table(as.character(args[1]),header=T)
s <- read.table(as.character(args[2]),header=T)

ind <- match(t$Rate_Category,s$Rate_Category)

t$Rate_Probability <- s$Rate_Probability[ind]
t$Rate_Factor <- s$Rate_Factor[ind]
t <- t %>% separate(Rate_Factor, into = c("dN", "dS"), sep = "\\,") 
t <- select(t,-c(Partition_Index, Rate_Category))

f=sub("pyvolve_out","final_site_rates",as.character(args[1]),fixed = TRUE)
write.table(t,f,sep="\t",quote=F,row.names=F)