library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

file_names = c("r4s_norm_rates","r4s_orig_rates")

for (name in file_names) {
  t1 <- list.files("mech_codon/r4s_rates/raw_rates",pattern=name,full.names=T)
  info = file.info(t1)
  t1 <- t1[info$size != 0]
  
  for (i in 1:length(t1)) 
  {
    r1 <- read.table(t1[i],skip=11,sep="\t")

    #reformat rate4site output
    r <- r1 %>% separate(V1,into=c("c1","c2","c3","c4","c5","c6","c7","c8"),sep="\\][:blank:]*|[:blank:]*\\[[:blank:]*|\\,[:blank:]*|[:blank:]+",extra="drop") %>%
      separate(c8,into=c("c8","c9"),sep="\\/") %>% 
      select(-c(c1,c9))
    colnames(r) <- c("pos","seq","score","qq_int_lower","qq_int_upper","stdev","num_taxa")
    
    str <- regexpr("_\\w+.txt",t1[i])[1]
    end <- regexpr(".txt",t1[i])[1]
    bias <- substr(t1[i],str+1,end-1)
    r$type <- rep(bias,length(r$num_taxa))
    
    str <- regexpr("bl\\d+",t1[i])[1]
    end <- regexpr("_\\w+.txt",t1[i])[1]
    bl <- as.numeric(substr(t1[i],str+2,end-1))
    r$bl <- rep(bl,length(r$num_taxa))
    
    str <- regexpr("rep\\d+",t1[i])[1]
    end <- regexpr("_n\\d+",t1[i])[1]
    rep <- as.numeric(substr(t1[i],str+3,end-1))
    r$rep <- rep(rep,length(r$num_taxa))
    
    if (name=="r4s_norm_rates") {
      f_type="r4s_norm"
    } else {
      f_type="r4s_orig"
    }
    true_rates_file_name <- gsub(paste0("r4s_rates/raw_rates/",f_type),"true_rates/sim",t1[i])
    t <- read.table(true_rates_file_name,header=T)
    ##get true dN/dS by solving for (ns_changes/ns_sites) / (s_changes/s_sites)
    t$dn.ds <- (t$ns_changes/t$ns_sites) / (t$s_changes/t$s_sites)
    r$true <- t$dn.ds
    r$score <- as.numeric(r$score)
    
    c <- cor.test(r$true,r$score,method = c("spearman"))

    r$cor <- rep(c$estimate,length(r$num_taxa))
    
    if (i==1) {
      d <- r
    } else d <- rbind(d, r)
  }  
  bias_r <- filter(d,type=="bias")
  nobias_r <- filter(d,type=="nobias")
  write.csv(bias_r,file=paste0("mech_codon/r4s_rates/processed_rates/all_",name,"_bias.csv"),quote=F)
  write.csv(nobias_r,file=paste0("mech_codon/r4s_rates/processed_rates/all_",name,"_nobias.csv"),quote=F)
}