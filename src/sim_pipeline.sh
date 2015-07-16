#!/bin/bash
#number of taxa
model=$1
#branch lengths
num_taxa=$2
#model for the simulation
br_len=$3
#number of the simulation
sim_num=$4

if [ ! -d "${model}" ]; then
	mkdir "${model}"
fi

#run the pipeline in the model's directory
cd ${model} 

tree=t${num_taxa}_b${br_len}_${sim_num}.tre ##tree file name
aln=seq_t${num_taxa}_b${br_len}_${sim_num}.fasta ##multiple sequence alignment file name
sim_rates=site_rates_t${num_taxa}_b${br_len}_${sim_num}.txt ##pyvolve output file name
sim_rates_info=site_rates_info_t${num_taxa}_b${br_len}_${sim_num}.txt ##pyvolve output file name
r4s_norm_rates=r4s_norm_rates_t${num_taxa}_b${br_len}_${sim_num}.txt ##rate4site output file name for norm rates

if [ ! -d "aln" ]; then
	mkdir "aln"
	mkdir "aln/nuc"
	mkdir "aln/aa"
fi

if [ ! -d "sim_site_rates/" ]; then
	mkdir "sim_site_rates/"
fi

##simulate multiple sequence alignment based on the tree
python ../src/simulate_aln.py $model trees/${tree} aln/nuc/$aln sim_site_rates/${sim_rates} sim_site_rates/${sim_rates_info}

##merge simulate_aln.py output
#Rscript ../src/merge_site_rates.r sim_site_rates/${sim_rates} sim_site_rates/${sim_rates_info}
# if [ -f "sim_site_rates/${sim_rates_info}" ]; then
# 	rm sim_site_rates/${sim_rates_info}
# fi

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