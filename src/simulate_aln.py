from pyvolve import *
import numpy as np
import sys

# Read in phylogeny
num = sys.argv[1]
bl = sys.argv[2]
tree = read_tree(file = sys.argv[3])
model_name=sys.argv[4].upper()
aln_file=sys.argv[5]

# Define model(s). State frequencies will be equal by default (1/61 per codon). You can tweak this by adding a vector, of length 61, w/ key "state_freqs" to the parameters dictionary. Should sum to 1 and is mapped to codons in alphabetical order: AAA, AAC, AAG, AAT, ...TTT
kappa = 4.5
parameters = {"kappa": kappa, "omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5

#gy_model = Model("GY", parameters, scale_matrix = "neutral") # add the rate_probs keyword argument **here** to specify different probabilities for dN/dS categories
mg_model = Model(model_name, parameters, scale_matrix = "neutral")


# Define partition(s)
length = 100 # number of codon positions
#gy_part = Partition(size = length, models = gy_model)
mg_part = Partition(size = length, models = mg_model)

# Evolve the GY and the MG sequences.
#gy_evolve = Evolver(partitions = gy_part, tree = tree)
#gy_evolve(seqfile = "../test_msa/nuc_msa/gy_sequences.fasta")

mg_evolve = Evolver(partitions = mg_part, tree = tree)
mg_evolve(seqfile = aln_file)
