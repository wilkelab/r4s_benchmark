# R code for generating a balanced tree
library(ape)

args <- commandArgs(trailingOnly = TRUE)
num <- as.numeric(args[1]) # num taxa must be a power of 2 for a binary tree ; 32, 64, 128, 256
tree <- stree(num, type = "balanced") # generated binary tree w/ 64 taxa. There are no branch lengths, though. Can add them by accessing tree$edge.length

bl <- as.numeric(args[2])
tree$edge.length <- rep(bl,nrow(tree$edge)) # sets all branch lengths to 0.01 ; 0.001, 0.0033, 0.01, 0.033, 0.1
##tree$edge.length<-runif(n=nrow(tree$edge),min=0,max=0.1) # draws branch lengths from uniform distribution with min0, max2

f = as.character(args[3])
write.tree(tree, file=f)
