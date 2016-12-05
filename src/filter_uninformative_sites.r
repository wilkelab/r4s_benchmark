library(tidyr)
library(dplyr)
library(bio3d)
library(readr)

setwd("r4s_benchmark/")

r <- read_csv("natural_prot/processed_rates/all_orig_rates.csv")
d <- data.frame()
for (name in unique(r$prot_name)) {
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
  bool_v[new_sites] <- rep(TRUE,length(new_sites))
  p$uninformative_site <- bool_v
  print(sites1)
  print(new_sites)
  write.csv(p,file=paste0("natural_prot/inferred_rates/",name,"_FEL1_reformatted.txt"))
}
