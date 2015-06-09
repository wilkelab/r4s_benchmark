from pyvolve import *
import numpy as np
import sys

# Read in phylogeny
model = sys.argv[1]
tree = read_tree(file = sys.argv[2])
aln_file=sys.argv[3]

# Define partition(s)
length = 100 # number of codon positions

if model=="dN" or model=="dN_dS": #aa mutation is symmetric 
	kappa = 4.5 #set transition:transversion ratio
	model_name="MG" #set model to mg-style
	
	if model=="dN":
		parameters = {"kappa": kappa, "omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5
	if model=="dN_dS":
		parameters = {"kappa": kappa, "alpha": np.arange(0.1, 1.6, 0.1), "beta": np.arange(0.1, 1.6, 0.1) }

	mg_model = Model(model_name, parameters, scale_matrix = "neutral")
	mg_part = Partition(size = length, models = mg_model)
	mg_evolve = Evolver(partitions = mg_part, tree = tree)
	mg_evolve(seqfile = aln_file)
	
#elif model="ms_dS" or model=="ms_no_dS":
#	model_name="MutSel"
else:
	sys.exit("wrong input model")
	


