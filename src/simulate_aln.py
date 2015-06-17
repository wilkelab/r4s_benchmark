from pyvolve import *
import numpy as np
import sys
from dnds_functions import *

# Read in phylogeny
model = sys.argv[1]
tree = read_tree(file = sys.argv[2])
aln_file=sys.argv[3]

# Define partition(s)
length = 100 # number of codon positions

if model=="dN" or model=="dN_dS": #aa mutation is symmetric 
	kappa = 4.5 #set transition:transversion ratio
	
	if model=="dN": ##varying dN 
		parameters = {"kappa": kappa, "omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5
	if model=="dN_dS": ##varying dN and dS
		##make every combination of rates 0.1-1.5 of dN to 0.1-1.5 of dS, i.e. 0.1 to 0.1, 0.1 to 0.2, etc.
		r = np.arange(0.1, 1.6, 0.1)
		n = len(r)
		l1 = np.repeat(r,n)
		l2 = np.tile(r,n)
		parameters = {"kappa": kappa, "alpha": l1, "beta": l2 }
	model = Model("MG", parameters, scale_matrix = "neutral")

# ######STEPHANIE'S CODE############
# from pyvolve import *
# import numpy as np
# 
# size = 100 # 100 sites to evolve
# t = read_tree("tre.txt")
# parts = []
# for i in range(size):
#     simulated_fitness = np.random.exponential(scale=1, size = 20) # draw 20 fitness values, so each codon family gets the same fitness. you should probably save these to a file.. the scale argument (i think that's the argument?) is the mean of the distribution. look this up, though..
#     m = Model("mutsel", {"fitness":simulated_fitness})     
#     p = Partition(models = m, size = 1)
#     parts.append(p)
#
#	  # dN/dS calculation
#dnds = None
#while dnds is None:
#     codon_freqs = dict(zip(codons, m.params["state_freqs"]))
#     mu = m.params["mu"]
#     try:
#		dnds = derive_dnds(codon_freqs, mu)
#	  except:
#		dnds = NA
# 
# e = Evolver(partitions = parts, tree = t)
# e()	
# ###################################

# elif model=="ms_dS" or model=="ms_no_dS":
# 	nuc_list = ["AC","CA","AG","GA","AT","TA","CG","GC","CT","TC","GT","TG"]
# 	rates =[ 1.33772085,  1.03590118,  0.1616799 ,  0.48613279,  1.43143935, 0.0910265 ,  0.56344429,  0.55420064,  0.46036002,  1.2304828 ,0.76040527,  0.30253117]
# 	mu = dict(zip(nuc_list,rates))
# 	#set all mut rates to 1
# 	
# 	codon_fitness = np.random.normal(size = 61) # constructs a vector of normally distributed codon fitness values, as an example
# 	##different fitness for each amino acid, but no dS
# 	##Stephanie uses exponential distr
# 	
# 	if model=="ms_dS":
# 		parameters = {"mu":mu,"fitness": codon_fitness}
#  	
# 	if model=="ms_no_dS":
#  		parameters = {"mu":mu,"fitness": codon_fitness}
# 		
# 	model = Model("MutSel",parameters)	

else:
	sys.exit("wrong input model")
	
part = Partition(size = length, models = model)
evolve = Evolver(partitions = part, tree = tree)
evolve(seqfile = aln_file)


