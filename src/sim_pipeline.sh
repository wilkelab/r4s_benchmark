#!/bin/bash
#number of taxa
num=$1
#branch lengths
bl=$2
tree=test_t${num}_b${bl}.tre
model=mg
aln=${model}_seq_t${num}_b${bl}.fasta

if [ ! -d "../trees" ]; then
	mkdir "../trees"
fi

##simulate a tree 
Rscript generate_balanced_tree_ape.R $num $bl ../trees/${tree}

if [ ! -d "../aln" ]; then
	mkdir "../aln"
	mkdir "../aln/nuc"
	mkdir "../aln/aa"
fi

##simulate multiple sequence alignment based on the tree
python simulate_aln.py $num $bl ../trees/${tree} $model ../aln/nuc/$aln

if [ ! -d "../sim_site_rates/" ]; then
	mkdir "../sim_site_rates/"
	mkdir "../sim_site_rates/pyvolve_out/"
	mkdir "../sim_site_rates/final_site_rates/"
fi

##move simulate_aln.py output
if [ -f "site_rates.txt" ]; then
	mv site_rates.txt ../sim_site_rates/pyvolve_out/site_rates_t${num}_b${bl}.txt
fi

if [ -f "site_rates_info.txt" ]; then
	mv site_rates_info.txt ../sim_site_rates/pyvolve_out/site_rates_info_t${num}_b${bl}.txt
fi

##merge pyvolve output
Rscript merge_site_rates.r ../sim_site_rates/pyvolve_out/site_rates_t${num}_b${bl}.txt ../sim_site_rates/pyvolve_out/site_rates_info_t${num}_b${bl}.txt 

##convert an alignment from nuc to aa 
python translate_aln.py ../aln/nuc/$aln

##run rate4site 
