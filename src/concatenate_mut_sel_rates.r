library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readr)

setwd("r4s_benchmark/")
file_names = c("r4s_norm_rates","r4s_orig_rates")
model="mut_sel"

t1_true_rates <- read_csv("../dnds_1rate_2rate/postprocessing/dataframes/substitution_counts.csv")
t_true_rates <- t1_true_rates %>% filter(pitype=="unequalpi")

inf_rates_nobias <- read_csv("../dnds_1rate_2rate/postprocessing/dataframes/results_balancedtrees_nobias_unequalpi.csv")
inf_rates_bias <- read_csv("../dnds_1rate_2rate/postprocessing/dataframes/results_balancedtrees_bias_unequalpi.csv")
fel1_nobias <- inf_rates_nobias %>% filter(method=="FEL1")
fel1_bias <- inf_rates_bias %>% filter(method=="FEL1")
fel1_r <- rbind(fel1_nobias,fel1_bias)

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
    biastype_str <- substr(t1[i],str+1,end-1)
    r$biastype <- rep(biastype_str,length(r$num_taxa))

    str <- regexpr("bl\\d+",t1[i])[1]
    end <- regexpr("_\\w+.txt",t1[i])[1]
    bl_num <- as.numeric(substr(t1[i],str+2,end-1))
    r$bl <- rep(bl_num,length(r$num_taxa))

    str <- regexpr("rep\\d+",t1[i])[1]
    end <- regexpr("_n\\d+",t1[i])[1]
    rep_num <- as.numeric(substr(t1[i],str+3,end-1))
    r$rep <- rep(rep_num,length(r$num_taxa))

    str <- regexpr("n\\d+",t1[i])[1]
    end <- regexpr("_bl\\d+",t1[i])[1]
    n <- as.numeric(substr(t1[i],str+1,end-1))
    
    num_taxa <- r$num_taxa[1]
  
    inf_r <- fel1_r %>% filter(rep==rep_num,bl==bl_num,ntaxa==num_taxa,biastype==biastype_str)
    r$inferred <- inf_r$dnds
    
    true_r <- t_true_rates %>% filter(rep==rep_num,bl==bl_num,ntaxa==num_taxa,biastype==biastype_str)
    r$true <- true_r$truednds
    
    if (name=="r4s_norm_rates") {
      f_type="r4s_norm"
    } else {
      f_type="r4s_orig"
    }

    if (i==1) {
      d <- r
    } else d <- rbind(d, r)
  }  
  #bias_r <- filter(d,biastype=="bias")
  #nobias_r <- filter(d,biastype=="nobias")
  #write.csv(bias_r,file=paste0(model,"/processed_rates/all_",name,"_bias.csv"),quote=F)
  #write.csv(nobias_r,file=paste0(model,"/processed_rates/all_",name,"_nobias.csv"),quote=F)
}