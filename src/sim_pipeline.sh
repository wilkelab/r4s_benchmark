#!/bin/bash
#number of taxa
model=$1
#branch lengths
sim_num=$2
#number of the simulation
num_taxa=$3
#model for the simulation
br_len=$4

if [ ! -d "${model}" ]; then
	mkdir "${model}"
fi

#run the pipeline in the model's directory
cd ${model}

tree=t${num_taxa}_b${br_len}_${sim_num}.tre ##tree file name
aln=${style}_seq_t${num_taxa}_b${br_len}.fasta ##multiple sequence alignment file name
sim_rates=site_rates_t${num_taxa}_b${br_len}_${sim_num}.txt ##pyvolve output file name
sim_rates_info=site_rates__info_t${num_taxa}_b${br_len}_${sim_num}.txt ##pyvolve output file name
r4s_norm_rates=r4s_norm_rates_t${num_taxa}_b${br_len}_${sim_num}.txt ##rate4site output file name for norm rates

if [ ! -d "trees" ]; then
	mkdir "trees"
fi

if [ -f "trees/${tree}" ]; then
	rm trees/${tree}
fi

##simulate a tree 
Rscript ../src/generate_balanced_tree_ape.R $num_taxa $br_len trees/${tree}

if [ ! -d "aln" ]; then
	mkdir "aln"
	mkdir "aln/nuc"
	mkdir "aln/aa"
fi

if [ ! -d "sim_site_rates/" ]; then
	mkdir "sim_site_rates/"
fi

if [ -f "aln/nuc/$aln" ]; then
	rm aln/nuc/$aln
fi

if [ -f "sim_site_rates/pyvolve_out/${sim_rates}" ]; then
	rm sim_site_rates/pyvolve_out/${sim_rates}
	rm sim_site_rates/pyvolve_out/${sim_rates_info}
	rm sim_site_rates/final_site_rates/${sim_rates}
fi

##simulate multiple sequence alignment based on the tree
python ../src/simulate_aln.py $model trees/${tree} aln/nuc/$aln

##move simulate_aln.py output
if [ -f "site_rates.txt" ]; then
	mv site_rates.txt sim_site_rates/${sim_rates}
	mv site_rates_info.txt sim_site_rates/${sim_rates_info}
fi

##merge simulate_aln.py output
Rscript ../src/merge_site_rates.r sim_site_rates/${sim_rates} sim_site_rates/${sim_rates_info}

if [ -f "sim_site_rates/${sim_rates_info}" ]; then
	rm sim_site_rates/${sim_rates_info}
fi

if [ -f "aln/aa/$aln" ]; then
	rm aln/aa/$aln
fi

##convert an alignment from nuc to aa 
python ../src/translate_aln.py aln/nuc/$aln

if [ ! -d "r4s_site_rates/" ]; then
	mkdir "r4s_site_rates/"
fi

##run rate4site 
../../rate4site.3.2.source/sourceMar09/rate4site -s aln/aa/$aln -t trees/$tree -o $r4s_norm_rates 
if [ -f r4s.res ]; then
	rm r4s.res
	rm r4sOrig.res
	rm TheTree.txt
	mv $r4s_norm_rates r4s_site_rates/
fi