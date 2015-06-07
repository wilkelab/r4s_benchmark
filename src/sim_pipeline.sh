#!/bin/bash
#number of taxa
num=$1
#branch lengths
bl=$2
tree=test_t${num}_b${bl}.tre
model=mg
aln=${model}_seq_t${num}_b${bl}.fasta
r4s_norm_rates=r4s_rates_t${num}_b${bl}.txt
r4s_unnorm_rates=r4s_rates_t${num}_b${bl}.txt
r4s_tree=r4s_rates_t${num}_b${bl}.txt

if [ -f "trees/${tree}" ]; then
	rm trees/${tree}
fi

if [ ! -d "trees" ]; then
	mkdir "trees"
fi

##simulate a tree 
Rscript src/generate_balanced_tree_ape.R $num $bl trees/${tree}

if [ ! -d "aln" ]; then
	mkdir "aln"
	mkdir "aln/nuc"
	mkdir "aln/aa"
fi

if [ -f "aln/nuc/$aln" ]; then
	rm aln/nuc/$aln
fi

if [ -f "sim_site_rates/pyvolve_out/site_rates_t${num}_b${bl}.txt" ]; then
	rm sim_site_rates/pyvolve_out/site_rates_t${num}_b${bl}.txt
fi

if [ -f "sim_site_rates/pyvolve_out/site_rates_info_t${num}_b${bl}.txt" ]; then
	rm sim_site_rates/pyvolve_out/site_rates_t${num}_b${bl}.txt
fi

if [ -f "sim_site_rates/final_site_rates/site_rates_t${num}_b${bl}.txt" ]; then
	rm sim_site_rates/final_site_rates/site_rates_t${num}_b${bl}.txt
fi

##simulate multiple sequence alignment based on the tree
python src/simulate_aln.py $num $bl trees/${tree} $model aln/nuc/$aln

if [ ! -d "sim_site_rates/" ]; then
	mkdir "sim_site_rates/"
	mkdir "sim_site_rates/pyvolve_out/"
	mkdir "sim_site_rates/final_site_rates/"
fi

##move simulate_aln.py output
if [ -f "site_rates.txt" ]; then
	mv site_rates.txt sim_site_rates/pyvolve_out/site_rates_t${num}_b${bl}.txt
fi

if [ -f "site_rates_info.txt" ]; then
	mv site_rates_info.txt sim_site_rates/pyvolve_out/site_rates_info_t${num}_b${bl}.txt
fi

##merge pyvolve output
Rscript src/merge_site_rates.r sim_site_rates/pyvolve_out/site_rates_t${num}_b${bl}.txt sim_site_rates/pyvolve_out/site_rates_info_t${num}_b${bl}.txt 

if [ -f "aln/aa/$aln" ]; then
	rm aln/aa/$aln
fi

##convert an alignment from nuc to aa 
python src/translate_aln.py aln/nuc/$aln

if [ ! -d "r4s_site_rates/" ]; then
	mkdir "r4s_site_rates/"
	mkdir "r4s_site_rates/r4s_out/"
	mkdir "r4s_site_rates/final_site_rates/"
fi

##run rate4site 
echo "../rate4site.3.2.source/sourceMar09/rate4site -s aln/aa/$aln -t trees/$tree -o $r4s_norm_rates -y $r4s_unnorm_rates -x $r4s_tree"
if [ -f $r4s_norm_rates ]; then
	mv $r4s_norm_rates r4s_site_rates/final_site_rates/
	mv $r4s_unnorm_rates r4s_site_rates/r4s_out/
	mv $r4s_tree r4s_site_rates/r4s_out/
fi