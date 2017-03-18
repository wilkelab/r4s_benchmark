# R code for generating a balanced tree
library(ape)

setwd("r4s_benchmark/")
num_taxa <- c(16, 32, 64, 128, 256, 512, 1024, 2048)
br_len <- c(0.0025,0.01,0.04,0.16,0.64)
for (n in num_taxa) {
  tree <- stree(n, type = "balanced") # generated binary tree w/ num of taxa.
  
  for (bl in br_len) {
    tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to bl
  
    n_2base=log(n,base=2)
    f = paste0("../r4s_benchmark_data/trees/n",n_2base,"_bl",bl,".tre")
    write.tree(tree, file=f)
  }
}
