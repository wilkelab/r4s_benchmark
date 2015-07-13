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

	part = Partition(size = length, models = model)
	evolve = Evolver(partitions = part, tree = tree)
	evolve(ratefile = rate_file, infofile = rate_info_file, seqfile = aln_file)

elif model=="ms_dS" or model=="ms_no_dS":
	parts = []
	
	base_name = rate_file.split("_")
	dnds_file_name = "dnds_site_rates_" + "_".join(base_name[-3:])
	dnds_file = open("sim_site_rates/"+dnds_file_name,"w")
	dnds_file.write("Site_Index\tdN/dS\n")
	dnds_list = []

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
	evolve(ratefile = rate_file, infofile = rate_info_file, seqfile = aln_file)	

else:
	sys.exit("wrong input model")
	



