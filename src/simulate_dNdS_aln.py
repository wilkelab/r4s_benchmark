from pyvolve import *
import numpy as np
import sys
import random
import math
import count_simulated_dnds
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def make_mc_model(bias, tree_file, aln_file, sim_rates_file, sim_rates_info_file, rate_distr):
	tree=read_tree(file = tree_file)
	kappa=4.5
	length=100
	
	if rate_distr=='uniform':
		dNdS_lst = np.random.uniform(0.1,1.6,10000)
		
	if rate_distr=='gamma':
		##shape and rate parameters used from A. G. Meyer, C. O. Wilke (2015). The utility of protein structure as a predictor of site-wise dN/dS varies widely among HIV-1 proteins.
		##gamma shape=0.312 and gamma rate=1.027 for integrase (Table 2 in the paper).
		##Set shape=0.312 and scale=1/rate=1/1.027. Mean=0.312*1/1.027= 0.304
		dNdS_lst = np.random.gamma(0.312,1/1.027,10000)
		
	if bias=="nobias": ##simulate varying dN 
		parameters = {"kappa":4.5, "omega": dNdS_lst} #set dN value range 0.1 - 1.6
	elif bias=="bias": ##simulate varying dN and dS
		##set dS value range 0.5 - 2
		if rate_distr=='uniform':
			dS_lst = np.random.uniform(0.5,2,10000)

		parameters = {"kappa":4.5, "alpha": dS_lst, "beta": dNdS_lst*dS_lst }
	else:
		print """
		Incorrect model name!
		"""
		
	model = Model("MG", parameters)
	parts = Partition(models = model, size = length)

	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = sim_rates_file, infofile = sim_rates_info_file, seqfile = aln_file, seqfmt = "fasta")
	
def main(argv):

	if len(argv) != 6: # wrong number of arguments
		print """Usage:
simulate_aln.py <bias> <tree_file> <aln_file> <sim_rates_file> <sim_rates_info_file> 
"""
		sys.exit()

	bias=argv[1]
	tree_file = argv[2]
	aln_file=argv[3]
	sim_rates_file=argv[4]
	sim_rates_info_file=argv[5]
	rate_distr=argv[6]
	
	make_mc_model(bias, tree_file, aln_file, sim_rates_file, sim_rates_info_file, rate_distr)
		
if __name__ == "__main__":
	main(sys.argv)