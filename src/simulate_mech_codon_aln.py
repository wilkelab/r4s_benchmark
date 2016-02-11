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

def make_mc_model(bias, tree_file, aln_file, ancestral_aln_file, sim_rates_file, sim_rates_info_file, true_rates_file):
	tree=read_tree(file = tree_file)
	kappa=4.5
	length=100

	if bias=="nobias": ##varying dN 
		parameters = {"kappa":4.5, "omega": np.arange(0.1, 1.6, 0.1) } # dN values ranging from 0.1 - 1.5
	elif bias=="bias": ##varying dN and dS
		##make every combination of rates 0.1-1.5 of dN to 0.1-1.5 of dS, i.e. 0.1 to 0.1, 0.1 to 0.2, etc.
		r = np.arange(0.1, 1.6, 0.1)
		n = len(r)
		l1 = np.repeat(r,n)
		l2 = np.tile(r,n)
		parameters = {"kappa":4.5, "alpha": l1, "beta": l2 }
	else:
		print """
		Incorrect model name!
		"""
		
	model = Model("MG", parameters)
	parts = Partition(models = model, size = length)

	evolve = Evolver(partitions = parts, tree = tree)
	evolve(ratefile = sim_rates_file, infofile = sim_rates_info_file, seqfile = aln_file)
	
	ancestral_aln_dict = evolve.get_sequences(anc = True) ##get simulated ancestor sequences 
	aln_list=[]
	
	for node in ancestral_aln_dict.keys():
		seq_r = SeqRecord(Seq(ancestral_aln_dict[node]), id=node, description='') #create a SeqRecord object from amino acid Seq object to start a list 
		aln_list.append(seq_r)
	msa = MultipleSeqAlignment(aln_list) # create an MSA object
	AlignIO.write(msa, ancestral_aln_file, "fasta") 
		
	# Mutation rate dictionary - FYI you can also grab this from pyvolve directly, just use the next line:
	mu = 1
	mu_dict = {'AT': mu, 'CG': mu, 'AC': mu, 'GT':mu, 'AG': kappa*mu,'GA': kappa*mu, 'TA': mu, 'GC': mu, 'CA': mu, 'TG':mu, 'GA': kappa*mu,'CT': kappa*mu, 'TC': kappa*mu}
	
	# Count dN/dS
	c = count_simulated_dnds.dNdS_Counter(ancestral_aln_file, tree_file, mu_dict)
	c.calculate_dnds(savefile = true_rates_file)

def main(argv):

	if len(argv) != 7: # wrong number of arguments
		print """Usage:
simulate_aln.py <bias> <tree_file> <aln_file> <sim_rates_file> <sim_rates_info_file> <ture_rates_file> 
"""
		sys.exit()

	bias=argv[1]
	tree_file = argv[2]
	aln_file=argv[3]
	sim_rates_file=argv[4]
	sim_rates_info_file=argv[5]
	true_rates_file=argv[6]
	
	ancestral_aln_file=aln_file.replace("nuc","ancestral")
	make_mc_model(bias,tree_file, aln_file, ancestral_aln_file, sim_rates_file, sim_rates_info_file, true_rates_file)
		
if __name__ == "__main__":
	main(sys.argv)