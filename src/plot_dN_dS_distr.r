library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
  t2 <- list.files(paste0(model,"/sim_site_rates/merged_output"),
                   full.names=T)
}

if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
  t2 <- list.files(paste0(model,"/sim_site_rates/inferred_rates"),
                   full.names=T)
}

for (i in 1:length(t2)) 
{
  if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
    s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","rate_probability","dN","dS"))
    s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS) 
  }
  if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
    s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","dN","dS"))
    s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS)
  }
  
  str <- regexpr("b\\d+",t2[i])[1]
  end <- regexpr("_\\d+.txt",t2[i])[1]
  bl <- as.numeric(substr(t2[i],str+1,end-1))
  print(bl)
  str <- regexpr("t\\d+",t2[i])[1]
  end <- regexpr("_b\\d+",t2[i])[1]
  nt <- as.numeric(substr(t2[i],str+1,end-1))
  
  str <- regexpr("_\\d+.txt",t2[i])[1]
  end <- regexpr(".txt",t2[i])[1]
  sim_num <- as.numeric(substr(t2[i],str+1,end-1))
  
  s$branch_len <- rep(bl,length(s$dN.dS))
  s$num_taxa <- rep(nt,length(s$dN.dS))
  s$sim_num <- rep(sim_num,length(s$dN.dS))
  
  if (i==1) {
    d <- s
  } else d <- rbind(d, s)
}

d <- d %>% arrange(num_taxa,branch_len)

p1 <- ggplot(d,aes(dN.dS)) + 
  geom_density()+
  facet_grid(num_taxa~branch_len)+
  xlab("dN/dS") 
ggsave(paste0("plots/",model,"_rate_distr.png"))


