library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

file_names = c("r4s_norm_rates","r4s_orig_rates")
model="mut_sel"

t_bias <- read.csv("../sitewise_dnds_mutsel/results/processed_results/full_results_bias_dataset.csv")
t_nobias <- read.csv("../sitewise_dnds_mutsel/results/processed_results/full_results_nobias_dataset.csv")
t_bias_filtered <- filter(t_bias,method=="FEL1")
t_nobias_filtered <- filter(t_nobias,method=="FEL1")

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
    
    if (bias=="bias") { 
      t_rate <- filter(t_bias_filtered,
                          bl==r$bl[1],
                          rep==r$rep[1],
                          ntaxa==r$num_taxa[1])
      r$true <- as.numeric(t_rate$true1)
    } else {
        t_rate <- filter(t_nobias_filtered,
                         bl==r$bl[1],
                         rep==r$rep[1],
                         ntaxa==r$num_taxa[1])
        r$true <- as.numeric(t_rate$true1)
    }
    
    r$score <- as.numeric(r$score)
    c <- cor.test(r$true,r$score,method = c("spearman"))
    
    r$cor_true <- rep(c$estimate,length(r$num_taxa))
    
    if (name=="r4s_orig_rates") {
      r$true_norm <- r$true/mean(r$true)
      r$score_norm <- r$score/mean(r$score)
      r$rmsd <- sqrt(sum((r$score_norm-r$true_norm)^2)/length(r$score_norm))
      r$bias <- r$score_norm-r$true_norm
    }
   
    if (i==1) {
      d <- r
    } else d <- rbind(d, r)
  }  
  bias_r <- filter(d,type=="bias")
  nobias_r <- filter(d,type=="nobias")
  write.csv(bias_r,file=paste0(model,"/r4s_rates/processed_rates/all_",name,"_bias.csv"),quote=F)
  write.csv(nobias_r,file=paste0(model,"/r4s_rates/processed_rates/all_",name,"_nobias.csv"),quote=F)
}