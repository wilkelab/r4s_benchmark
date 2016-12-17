library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(bio3d)

file_names = c("norm_rates","orig_rates")

setwd("r4s_benchmark/")
for (name in file_names) {
  t1 <- list.files("natural_prot/r4s_rates/raw_rates",pattern=name)

  for (i in 1:length(t1)) 
  {
    r1 <- read.table(paste0("natural_prot/r4s_rates/raw_rates/",t1[i]),skip=11,sep="\t")
    #reformat rate4site output
    r <- r1 %>% separate(V1,into=c("c1","c2","c3","c4","c5","c6","c7","c8"),sep="\\][:blank:]*|[:blank:]*\\[[:blank:]*|\\,[:blank:]*|[:blank:]+",extra="drop") %>%
      separate(c8,into=c("c8","c9"),sep="\\/") %>% 
      select(-c(c1,c9))
    colnames(r) <- c("pos","seq","score","qq_int_lower","qq_int_upper","stdev","num_taxa")
    
    str <- 1
    end <- regexpr(paste0("_",name,".txt"),t1[i])[1]
    protein_name <- substr(t1[i],str,end-1)
    r$prot_name <- rep(protein_name,length(r$num_taxa))
    
    if (name=="r4s_norm_rates") {
      f_type="r4s_norm"
    } else {
      f_type="r4s_orig"
    }
    
    print(protein_name)
    if (substr(protein_name,1,4)=="ENST") {
      aln_file <- paste0("natural_prot/aln/processed_aln/",protein_name,"_p_aligned.fasta")
    } else {
      aln_file <- paste0("natural_prot/aln/processed_aln/",protein_name,"_clean_protein.fasta")
    }
    f <- read.fasta(aln_file)
    
    inferred_rates_file_name <- paste0("natural_prot/inferred_rates/",protein_name,"_FEL1_reformatted.txt")
    inf_r <- read.csv(inferred_rates_file_name)

    unchanged_sites_file_name <- paste0("natural_prot/filtered_sites/",protein_name,"_unchanged_sites.txt")
    sites_t <- read.table(unchanged_sites_file_name,header=T)
    
    filtered_inferred <- inf_r$dN.dS
    filtered_inferred[sites_t$unchanged_site]=rep(0,length(which(sites_t$unchanged_site)))

    r$inferred <- filtered_inferred[which(!is.gap(f$ali[1,]))]

    if (i==1) {
      d <- r
    } else d <- rbind(d, r)
  }  
  write.csv(d,file=paste0("natural_prot/processed_rates/all_",name,".csv"),quote=F)
}