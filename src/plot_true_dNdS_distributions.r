library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(readr)

setwd("r4s_benchmark/")
r_bias_dNdS_raw <- read_csv("mech_codon/processed_rates/all_r4s_orig_rates_bias.csv")
r_nobias_dNdS_raw <- read_csv("mech_codon/processed_rates/all_r4s_orig_rates_nobias.csv")
r_bias_MutSel_raw <- read_csv("mut_sel/processed_rates/all_r4s_orig_rates_bias.csv")
r_nobias_MutSel_raw <- read_csv("mut_sel/processed_rates/all_r4s_orig_rates_nobias.csv")

p_nobias_dNdS <- ggplot(r_nobias_dNdS_raw,aes(x=true,group=rep)) + 
    geom_density(alpha=0.3) +
    facet_grid(num_taxa~bl)

p_bias_dNdS <- ggplot(r_bias_dNdS_raw,aes(x=true,group=rep)) + 
  geom_density(alpha=0.3) +
  facet_grid(num_taxa~bl)

p_nobias_MutSel <- ggplot(r_nobias_MutSel_raw,aes(x=true,group=rep)) + 
  geom_density(alpha=0.3) +
  facet_grid(num_taxa~bl)

p_bias_MutSel <- ggplot(r_bias_MutSel_raw,aes(x=true,group=rep)) + 
  geom_density(alpha=0.3) +
  facet_grid(num_taxa~bl)
