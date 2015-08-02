from pyvolve import *
import numpy as np
import sys
from dnds_functions import *

# Read in phylogeny
model = sys.argv[1]
tree = read_tree(file = sys.argv[2])
aln_file=sys.argv[3]
rate_file=sys.argv[4]
rate_info_file=sys.argv[5]

# Define partition(s)
length = 100 # number of codon positions
parts=[]

if model=="dN" or model=="dN_dS": #aa mutation is symmetric 
	# mu={"AC'":1,"AG":4.5,"AT":1,
# 	"CA":1,"CG":1,"CT":4.5,
# 	"GA":4.5,"GC":1,"GT":1,
# 	"TA":1,"TC":4.5,"TG":1} # is the same as kappa = 4.5 which sets transition:transversion ratio
	kappa=4.5
	
	if model=="dN": ##varying dN 
		parameters = {"kappa":4.5, "omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5
	if model=="dN_dS": ##varying dN and dS
		##make every combination of rates 0.1-1.5 of dN to 0.1-1.5 of dS, i.e. 0.1 to 0.1, 0.1 to 0.2, etc.
		r = np.arange(0.1, 1.6, 0.1)
		n = len(r)
		l1 = np.repeat(r,n)
		l2 = np.tile(r,n)
		parameters = {"kappa":4.5, "alpha": l1, "beta": l2 }
	
	model = Model("MG", parameters, scale_matrix = "neutral")
	parts = Partition(models = model, size = length)
	
	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = rate_file, infofile = rate_info_file, seqfile = aln_file)

elif model=="ms_dS" or model=="ms_no_dS":
	
	out_file=rate_file.replace('assigned_rates','inferred_rates')
	dnds_file=open(out_file,"w")
	dnds_file.write("Site_Index\tdN_dS\n")
	
	if model=="ms_dS": 
		for i in range(length):
			simulated_fitness = np.random.exponential(scale=1, size = 61) # draw 61 fitness values for each codon, scale is the mean or beta of the exponential distribution
			model = Model("mutsel", {"fitness":simulated_fitness})     
			p = Partition(models = model, size = 1)
			parts.append(p)
			
			##calculate dN/dS
			codon_freqs = dict(zip(codons, model.params["state_freqs"]))
			mu = model.params["mu"]
			try:
				dnds = derive_dnds(codon_freqs, mu)
			except:
				dnds = NA
			dnds_file.write("%d\t%f\n" % (i+1, dnds))
			
	if model=="ms_no_dS":
		for i in range(length):
			 simulated_fitness = np.random.exponential(scale=1, size = 20) # draw 20 fitness values for amino acid, scale is the mean or beta of the exponential distribution
			 model = Model("mutsel", {"fitness":simulated_fitness})     
			 p = Partition(models = model, size = 1)
			 parts.append(p)
			 
			 codon_freqs = dict(zip(codons, model.params["state_freqs"]))
			 mu = model.params["mu"]
			 try:
			 	dnds = derive_dnds(codon_freqs, mu)
			 except:
			 	dnds = NA
			 dnds_file.write("%d\t%f\n" % (i+1, dnds))

	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = None, infofile = None, seqfile = aln_file)	

else:
	sys.exit("wrong input model")
	



