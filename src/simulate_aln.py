from pyvolve import *
import numpy as np
import sys
import dnds_functions
import random
import math
import slaculator 

def make_mc_model(var, tree_file, aln_file, rate_file, rate_info_file, length, inf_rate_file):
	tree=read_tree(file = tree_file)
	kappa=4.5

	if var=="dN": ##varying dN 
		parameters = {"kappa":4.5, "omega": np.arange(0.1, 1.6, 0.1) } # dN/dS values ranging from 0.1 - 1.5
	if var=="dN_dS": ##varying dN and dS
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
	
	# Mutation rate dictionary - FYI you can also grab this from pyvolve directly, just use the next line:
	# mu = 1
	mu_dict = {'AT': mu, 'CG': mu, 'AC': mu, 'GT':mu, 'AG': kappa*mu,'GA': kappa*mu, 'TA': mu, 'GC': mu, 'CA': mu, 'TG':mu, 'GA': kappa*mu,'CT': kappa*mu, 'TC': kappa*mu}
	
	# Count dN/dS
	print "Counting site-specific dN/dS"
	my_slaculator = slaculator.Slaculator(aln_file, tree_file, mu_dict)
	my_slaculator.calculate_dnds() #savefile = inf_rate_file) # To change name of output file from default 'slaculator_output.txt', add argument savefile = "filename.txt". It will be tab-delimited.

def make_ms_model(var, tree_file, aln_file, rate_file, rate_info_file, length, inf_rate_file):
	tree=read_tree(file = tree_file)
	dnds_file=open(inf_rate_file,"w")
	dnds_file.write("Site_Index\tdN\tdS\n")

	parts = []			
	if var=="dN":
		for i in range(length):
			simulated_fitness = np.random.exponential(scale=1, size = 20) # draw 20 fitness values for amino acid, scale is the mean or beta of the exponential distribution
			model = Model("mutsel", {"fitness":simulated_fitness})     
			p = Partition(models = model, size = 1)
			parts.append(p)

			codon_freqs = dict(zip(dnds_functions.codons, model.params["state_freqs"]))
			mu = model.params["mu"]
			
			try:
				dn,ds = dnds_functions.derive_dnds(codon_freqs, mu)
			except:
				dn = None
				ds = None
			dnds_file.write("%d\t%s\t%s\n" % (i+1, dn, ds))
				
	if var=="dN_dS": 	
		codon_index = {'A': [37, 38, 39, 40], 'C': [55, 57], 'E': [33, 35], 'D': [34, 36], 'G': [41, 42, 43, 44], 'F': [59, 61], 'I': [13, 14, 16], 'H': [18, 20], 'K': [1, 3], 'M': [15], 'L': [29, 30, 31, 32, 58, 60], 'N': [2, 4], 'Q': [17, 19], 'P': [21, 22, 23, 24], 'S': [10, 12, 51, 52, 53, 54], 'R': [9, 11, 25, 26, 27, 28], 'T': [5, 6, 7, 8], 'W': [56], 'V': [45, 46, 47, 48], 'Y': [49, 50]}
		index_list = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]

		mu = 1.
		kappa = 1.
		mu_dict = {'AT': mu, 'CG': mu, 'AC': mu, 'GT':mu, 'AG': kappa*mu,'GA': kappa*mu, 'TA': mu, 'GC': mu, 'CA': mu, 'TG':mu, 'GA': kappa*mu,'CT': kappa*mu, 'TC': kappa*mu}

		for i in range(length):
			'''
			step 1
			aa freqs from lambda distribution
			1. randomly draw lambda from some distribution
			2. randomly shuffle aa index
			'''
			lamb = np.random.uniform(0.2, 1.2)  #lambda value
			raw_freq_aa = [math.exp(-x*lamb) for x in index_list] #raw aa freqs list
			norm_freq_aa = [x/sum(raw_freq_aa) for x in raw_freq_aa] #normalize to 1
			random.shuffle(norm_freq_aa) #shuffle the freq list
			aa_freq =  dict(zip(codon_index.keys(), norm_freq_aa))
			'''
			step 2
			assign codon bias to synonymous codons
			''' 
			codon_freq = [0]*61
			codon_bias = {x:0 for x in codon_index.keys()} 
			for aa in codon_index.keys():
				num = len(codon_index[aa]) #number of synonymous codons
				lamb_cod = np.random.uniform(0.1,3) ##assigning weak bias
				raw_freq_cod = [math.exp(-x*lamb_cod) for x in range(num)]# draw #synonymous_codons distribution from exponential as raw 
				norm_freq_cod = [aa_freq[aa]*x/sum(raw_freq_cod) for x in raw_freq_cod] #normalize to aa freq
				for j in range(num):
					codon_freq[codon_index[aa][j]-1] = norm_freq_cod[j] #for codon index starts from 1
				codon_bias[aa] = np.var(norm_freq_cod)/np.mean(norm_freq_cod)  # codon bias in coefficient of variance
			
			freq_dict = dict(zip(dnds_functions.codons, codon_freq))

			model = Model("mutsel", {"state_freqs": codon_freq, "mu":mu_dict})
			p = Partition(models = model, size = 1)
			parts.append(p)
		 
			mu = model.params["mu"]			
			try:
				dn,ds = dnds_functions.derive_dnds(freq_dict, mu)
			except:
				dn = None
				ds = None
			dnds_file.write("%d\t%f\t%f\n" % (i+1, dn, ds))
				
	##Evolving sequences
	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = None, infofile = None, seqfile = aln_file)

def main(argv):

	if len(argv) != 6: # wrong number of arguments
		print """Usage:
simulate_aln.py <model> <tree_file> <aln_file> <rate_file> <rate_info_file> 
"""
		sys.exit()
		
	model = argv[1]
	tree_file = argv[2]
	aln_file=argv[3]
	rate_file=argv[4]
	rate_info_file=argv[5]

	# Define partition(s)
	length = 100 # number of codon positions
	
	##inferred rates file	
	inf_rate_file=rate_file.replace('assigned_rates','inferred_rates')
	
	if model=="mech_codon_dN":
		make_mc_model("dN",tree_file, aln_file, rate_file, rate_info_file, length, inf_rate_file)
	elif model=="mech_codon_dN_dS":
		make_mc_model("dN_dS",tree_file, aln_file, rate_file, rate_info_file, length, inf_rate_file)
	elif model=="mut_sel_dN":
		make_ms_model("dN",tree_file, aln_file, rate_file, rate_info_file, length, inf_rate_file)
	elif model=="mut_sel_dN_dS":
		make_ms_model("dN_dS",tree_file, aln_file, rate_file, rate_info_file, length, inf_rate_file)
	else:
		print """
		Incorrect model name!
		"""
		
if __name__ == "__main__":
	main(sys.argv)