library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t_bias <- read.table("codon_freq_lib_bias.txt",sep=" ")
t_nobias <- read.table("codon_freq_lib_nobias.txt",sep=" ")

codons=c("AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", 
         "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", 
         "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", 
         "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", 
         "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", 
         "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", 
         "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAC", 
         "TAT", "TCA", "TCC", "TCG", "TCT", "TGC", "TGG", 
         "TGT", "TTA", "TTC", "TTG", "TTT")
colnames(t_bias) <- codons
colnames(t_nobias) <- codons
fr_bias <- t_bias %>% gather(codons,frequency)
fr_nobias <- t_nobias %>% gather(codons,frequency)



