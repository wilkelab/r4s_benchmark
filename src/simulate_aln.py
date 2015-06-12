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
	
	if model=="dN": ##varying dN 
		parameters = {"kappa": kappa, "omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5
	if model=="dN_dS": ##varying dN and dS
		##make every combination of rates 0.1-1.5 of dN to 0.1-1.5 of dS, i.e. 0.1 to 0.1, 0.1 to 0.2, etc.
		r = np.arange(0.1, 1.6, 0.1)
		n = len(r)
		l1 = np.repeat(r,n)
		l2 = np.tile(r,n)
		parameters = {"kappa": kappa, "alpha": l1, "beta": l2 }

	mg_model = Model(model_name, parameters, scale_matrix = "neutral")
	mg_part = Partition(size = length, models = mg_model)
	mg_evolve = Evolver(partitions = mg_part, tree = tree)
	mg_evolve(seqfile = aln_file)
	
#elif model="ms_dS" or model=="ms_no_dS":
#	model_name="MutSel"
else:
	sys.exit("wrong input model")
	


