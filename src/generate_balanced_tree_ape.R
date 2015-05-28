# R code for generating a balanced tree
library(ape)
n <- 64 # num taxa must be a power of 2 for a binary tree
tree <- stree(n, type = "balanced") # generated binary tree w/ 64 taxa. There are no branch lengths, though. Can add them by accessing tree$edge.length
tree$edge.length <- 0.5 # sets all branch lengths to 0.5
tree$edge.length<-runif(n=**nrow(tree$edge),min=0,max=2) # draws branch lengths from uniform distribution with min0, max2
