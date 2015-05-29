from pyvolve import *
import numpy as np

# Read in phylogeny
tree = read_tree(file = "../test_trees/test1.tre")

# Define model(s). State frequencies will be equal by default (1/61 per codon). You can tweak this by adding a vector, of length 61, w/ key "state_freqs" to the parameters dictionary. Should sum to 1 and is mapped to codons in alphabetical order: AAA, AAC, AAG, AAT, ...TTT
kappa = 4.5
parameters = {"kappa": kappa, "omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5
#gy_model = CodonModel("GY", parameters, "neutral")
#gy_model.construct_model() # by default, each dN/dS category is equally likely. To change the probability of each dN/dS value occurring, add here the argument rate_probs = [list of probabilities associated with dN/dS categories specified in the parameters dictionary

mg_model = CodonModel("MG", parameters, "neutral")
mg_model.construct_model()

# Define partition(s)
length = 100 # number of codon positions
#gy_part = Partition(size = length, models = gy_model)
mg_part = Partition(size = length, models = mg_model)

# Evolve the GY and the MG sequences.
#gy_evolve = Evolver(partitions = gy_part, tree = tree)
#gy_evolve(seqfile = "gy_sequences.fasta")

mg_evolve = Evolver(partitions = mg_part, tree = tree)
mg_evolve(seqfile = "../test_msa/mg_sequences.fasta")
