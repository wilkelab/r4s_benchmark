library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)
model <- "mut_sel_dN"
t1 <- list.files(paste0(model,"/r4s_site_rates/orig_rates"),full.names=T)
info = file.info(t1)
t1 <- t1[info$size != 0]

if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
  t2 <- list.files(paste0(model,"/sim_site_rates/merged_output"),
                   full.names=T)
}

if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
  t2 <- list.files(paste0(model,"/sim_site_rates/inferred_rates"),
                   full.names=T)
}
t2 <- t2[info$size != 0]

for (i in 1:length(t1)) 
{
  print(t1[i])
  print(t2[i])
  
  r1 <- read.table(t1[i],skip=11,sep="\t")
  
  #reformat rate4site output
  r <- r1 %>% separate(V1,into=c("c1","c2","c3","c4","c5","c6","c7","c8"),sep="\\][:blank:]*|[:blank:]*\\[[:blank:]*|\\,[:blank:]*|[:blank:]+",extra="drop") %>%
    separate(c8,into=c("c8","c9"),sep="\\/") %>% 
    select(-c(c1,c9))
  colnames(r) <- c("pos","seq","score","qq_int_lower","qq_int_upper","stdev","num_taxa")
  
  if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
    s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","rate_probability","dN","dS"))
    s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS) 
  }
  if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
    s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","dN","dS"))
    s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS)
  }
  
  s$r4s_score <- as.numeric(r$score)
  s$num_taxa <- as.numeric(r$num_taxa)
  
  str <- regexpr("b\\d+",t1[i])[1]
  end <- regexpr("_\\d+.txt",t1[i])[1]
  bl <- as.numeric(substr(t1[i],str+1,end-1))
  print(bl)
  str <- regexpr("_\\d+.txt",t1[i])[1]
  end <- regexpr(".txt",t1[i])[1]
  sim_num <- as.numeric(substr(t1[i],str+1,end-1))
  
  s$branch_len <- rep(bl,length(s$num_taxa))
  s$sim_num <- rep(sim_num,length(s$num_taxa))
  
  if (model=="mut_sel_dN" | model=="mut_sel_dN_dS" | model=="mech_codon_dN_dS"){
    cor_test <- cor.test(s$dN.dS,s$r4s_score,method="spearman")
    cor <- cor_test$estimate
    p_val <- cor_test$p.value
  }
  
  if (model=="mech_codon_dN"){
    cor_test <- cor.test(s$dN,s$r4s_score,method="spearman")
    cor <- cor_test$estimate
    p_val <- cor_test$p.value
  }
  
  s$cor <- rep(cor,length(s$num_taxa))
  s$p_val <- rep(p_val,length(s$num_taxa))
  
  if (i==1) {
    d <- s
  } else d <- rbind(d, s)
}

d1 <- d %>% arrange(num_taxa,branch_len)

args <- commandArgs(trailingOnly = TRUE)
model <- "mut_sel_dN_dS"
t1 <- list.files(paste0(model,"/r4s_site_rates/orig_rates"),full.names=T)
info = file.info(t1)
t1 <- t1[info$size != 0]

if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
  t2 <- list.files(paste0(model,"/sim_site_rates/merged_output"),
                   full.names=T)
}

if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
  t2 <- list.files(paste0(model,"/sim_site_rates/inferred_rates"),
                   full.names=T)
}
t2 <- t2[info$size != 0]

for (i in 1:length(t1)) 
{
  print(t1[i])
  print(t2[i])
  
  r1 <- read.table(t1[i],skip=11,sep="\t")
  
  #reformat rate4site output
  r <- r1 %>% separate(V1,into=c("c1","c2","c3","c4","c5","c6","c7","c8"),sep="\\][:blank:]*|[:blank:]*\\[[:blank:]*|\\,[:blank:]*|[:blank:]+",extra="drop") %>%
    separate(c8,into=c("c8","c9"),sep="\\/") %>% 
    select(-c(c1,c9))
  colnames(r) <- c("pos","seq","score","qq_int_lower","qq_int_upper","stdev","num_taxa")
  
  if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
    s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","rate_probability","dN","dS"))
    s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS) 
  }
  if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
    s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","dN","dS"))
    s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS)
  }
  
  s$r4s_score <- as.numeric(r$score)
  s$num_taxa <- as.numeric(r$num_taxa)
  
  str <- regexpr("b\\d+",t1[i])[1]
  end <- regexpr("_\\d+.txt",t1[i])[1]
  bl <- as.numeric(substr(t1[i],str+1,end-1))
  print(bl)
  str <- regexpr("_\\d+.txt",t1[i])[1]
  end <- regexpr(".txt",t1[i])[1]
  sim_num <- as.numeric(substr(t1[i],str+1,end-1))
  
  s$branch_len <- rep(bl,length(s$num_taxa))
  s$sim_num <- rep(sim_num,length(s$num_taxa))
  
  if (model=="mut_sel_dN" | model=="mut_sel_dN_dS" | model=="mech_codon_dN_dS"){
    cor_test <- cor.test(s$dN.dS,s$r4s_score,method="spearman")
    cor <- cor_test$estimate
    p_val <- cor_test$p.value
  }
  
  if (model=="mech_codon_dN"){
    cor_test <- cor.test(s$dN,s$r4s_score,method="spearman")
    cor <- cor_test$estimate
    p_val <- cor_test$p.value
  }
  
  s$cor <- rep(cor,length(s$num_taxa))
  s$p_val <- rep(p_val,length(s$num_taxa))
  
  if (i==1) {
    d <- s
  } else d <- rbind(d, s)
}

d2 <- d %>% arrange(num_taxa,branch_len)
d1$model <- rep("mut_sel_dN",length(d1$pos))
d2$model <- rep("mut_sel_dN_dS",length(d2$pos))
d <- bind_rows(d1,d2)

p1 <- ggplot(d,aes(r4s_score),group=model) + 
  geom_density(aes(fill=model),alpha=0.4)+
  scale_x_continuous(breaks=seq(0,4,1), limits = c(0,4)) + 
  facet_grid(num_taxa~branch_len)+
  xlab("Rate4Site score") 
ggsave(paste0("plots/",model,"_r4s_score_distr.png"))

