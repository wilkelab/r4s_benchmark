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

def make_mc_model(bias, tree_file, aln_file, sim_rates_file, sim_rates_info_file):
	tree=read_tree(file = tree_file)
	kappa=4.5
	length=100
	
	dNdS_lst = np.random.uniform(0.1,1.6,10000)
	if bias=="nobias": ##simulate varying dN 
		parameters = {"kappa":4.5, "omega": dNdS_lst} #set dN value range 0.1 - 1.6
	elif bias=="bias": ##simulate varying dN and dS
		##set dS value range 0.5 - 2
		dS_lst = np.random.uniform(0.5,2,10000)
		parameters = {"kappa":4.5, "alpha": dS_lst, "beta": dNdS_lst*dS_lst }
	else:
		print """
		Incorrect model name!
		"""
		
	model = Model("MG", parameters)
	parts = Partition(models = model, size = length)

	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = sim_rates_file, infofile = sim_rates_info_file, seqfile = aln_file)
	
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
	
	make_mc_model(bias, tree_file, aln_file, sim_rates_file, sim_rates_info_file)
		
if __name__ == "__main__":
	main(sys.argv)