library(tidyr)
library(dplyr)
library(bio3d)
library(readr)

setwd("r4s_benchmark/")

r_nat_prot <- read_csv("natural_prot/processed_rates/all_orig_rates.csv")

d <- data.frame()
for (name in unique(r_nat_prot$prot_name)) {
  p <- read.csv(paste0("natural_prot/inferred_rates/",name,"_FEL1.txt"))
  sites1 <- which(p$dN.dS==1)
  
  if (grepl("hiv",name)) {
    aln_file <- paste0("natural_prot/aln/processed_aln/",name,"_clean_protein.fasta")
  } else aln_file <- paste0("natural_prot/aln/processed_aln/",name,"_p_aligned.fasta")

  aln_obj <- read.fasta(file = aln_file)
  aln <- aln_obj$ali
  rownames(aln) <- NULL
  
  new_sites <- c()
  for (i in sites1) {
    if (length(table(aln[,i]))==1) {
      new_sites <- c(new_sites,i)
    } else next
  }
  
  bool_v <- rep(FALSE,length(p$dN.dS))
  bool_v[new_sites] <- rep(0,length(new_sites))
  p$uninformative_site <- bool_v
  print(sites1)
  print(new_sites)
  write.csv(p,file=paste0("natural_prot/inferred_rates/",name,"_FEL1_reformatted.txt"))
}

r_mutsel_bias <- read_csv("mut_sel/processed_rates/all_r4s_norm_rates_bias.csv")
r_mutsel_nobias <- read_csv("mut_sel/processed_rates/all_r4s_norm_rates_nobias.csv")
r_dNdS_bias <- read_csv("mech_codon/processed_rates/all_r4s_orig_rates_bias.csv")
r_dNdS_nobias <- read_csv("mech_codon/processed_rates/all_r4s_orig_rates_nobias.csv")

d <- data.frame()
for (name in unique(r_nat_prot$prot_name)) {
  p <- read.csv(paste0("natural_prot/inferred_rates/",name,"_FEL1.txt"))
  sites1 <- which(p$dN.dS==1)

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
  
  str <- regexpr("n\\d+",t1[i])[1]
  end <- regexpr("_bl\\d+",t1[i])[1]
  n <- as.numeric(substr(t1[i],str+1,end-1))
  if (grepl("hiv",name)) {
    aln_file <- paste0("natural_prot/aln/processed_aln/",name,"_clean_protein.fasta")
  } else aln_file <- paste0("natural_prot/aln/processed_aln/",name,"_p_aligned.fasta")
  
  aln_obj <- read.fasta(file = aln_file)
  aln <- aln_obj$ali
  rownames(aln) <- NULL
  
  new_sites <- c()
  for (i in sites1) {
    if (length(table(aln[,i]))==1) {
      new_sites <- c(new_sites,i)
    } else next
  }
  
  bool_v <- rep(FALSE,length(p$dN.dS))
  bool_v[new_sites] <- rep(0,length(new_sites))
  p$uninformative_site <- bool_v
  print(sites1)
  print(new_sites)
  write.csv(p,file=paste0("natural_prot/inferred_rates/",name,"_FEL1_reformatted.txt"))
}
