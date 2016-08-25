library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)

file_names = c("r4s_norm_rates","r4s_orig_rates")
model="mut_sel"
t1_fel1 <- read_csv("../dnds_1rate_2rate/postprocessing/dataframes/full_results_equalpi_bias.csv") 
t2_fel1 <- read_csv("../dnds_1rate_2rate/postprocessing/dataframes/full_results_equalpi_nobias.csv") 
fel1_bias <- filter(t1_fel1,method=="FEL1")
fel1_nobias <- filter(t2_fel1,method=="FEL1")
fel1_t <- rbind(fel1_bias,fel1_nobias)

setwd("r4s_benchmark/")
for (name in file_names) {
  t1 <- list.files(paste0(model,"/r4s_rates/raw_rates"),pattern=name,full.names=T)
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
    type1 <- substr(t1[i],str+1,end-1)
    r$type <- rep(type1,length(r$num_taxa))
    
    str <- regexpr("bl\\d+",t1[i])[1]
    end <- regexpr("_\\w+.txt",t1[i])[1]
    bl1 <- as.numeric(substr(t1[i],str+2,end-1))
    r$bl <- rep(bl1,length(r$num_taxa))
    
    str <- regexpr("rep\\d+",t1[i])[1]
    end <- regexpr("_n\\d+",t1[i])[1]
    rep1 <- as.numeric(substr(t1[i],str+3,end-1))
    r$rep <- rep(rep1,length(r$num_taxa))
    
    str <- regexpr("n\\d+",t1[i])[1]
    end <- regexpr("_bl\\d+",t1[i])[1]
    num1 <- as.numeric(substr(t1[i],str+1,end-1))
    
    num_taxa1 <- r$num_taxa[1]
   
    fel1_r <- filter(fel1_t,rep==rep1,bl==bl1,ntaxa==num_taxa1,type==type1)
    r$inferred <- fel1_r$dnds
    r$true <- fel1_r$true
    
    if (name=="r4s_norm_rates") {
      f_type="r4s_norm"
    } else {
      f_type="r4s_orig"
    }
    
    r$score <- as.numeric(r$score)
    
    if (name=="r4s_orig_rates") {
      r$true_norm <- r$true/mean(na.omit(r$true))
      r$score_norm <- r$score/mean(r$score)
      r$inferred_norm <- r$inferred/mean(na.omit(r$inferred))
    }
   
    if (i==1) {
      d <- r
    } else d <- rbind(d, r)
  }  
  bias_r <- filter(d,type=="bias")
  nobias_r <- filter(d,type=="nobias")
  write.csv(bias_r,file=paste0(model,"/processed_rates/all_",name,"_bias.csv"),quote=F)
  write.csv(nobias_r,file=paste0(model,"/processed_rates/all_",name,"_nobias.csv"),quote=F)
}