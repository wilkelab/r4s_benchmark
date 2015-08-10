#!/bin/bash
#number of taxa
model=$1
#branch lengths
num_taxa=$2
#model for the simulation
br_len=$3
#number of the simulation
sim_num=$4


aln=seq_t${num_taxa}_b${br_len}_${sim_num}.fasta ##multiple sequence alignment file name
r4s_norm_rates=r4s_norm_rates_t${num_taxa}_b${br_len}_${sim_num}.txt ##rate4site output file name for norm rates
tree=t${num_taxa}_b${br_len}.tre

if [ ! -d "$SCRATCH/r4s_benchmark_data/${model}/aln/aa/" ]; then
	mkdir "$SCRATCH/r4s_benchmark_data/${model}/aln/aa/"
fi

##convert an alignment from nuc to aa 
python $HOME/r4s_benchmark/src/translate_aln.py $SCRATCH/r4s_benchmark_data/${model}/aln/nuc/$aln $SCRATCH/r4s_benchmark_data/${model}/aln/aa/$aln

if [ ! -d "$HOME/r4s_benchmark/${model}/r4s_site_rates/" ]; then
	mkdir "$HOME/r4s_benchmark/${model}/r4s_site_rates/"
fi

##run rate4site 
rate4site.3.2.source/sourceMar09/rate4site -s $SCRATCH/r4s_benchmark_data/${model}/aln/aa/${aln} -t $SCRATCH/r4s_benchmark_data/trees/${tree} -o $HOME/r4s_benchmark/${model}/r4s_site_rates/${r4s_norm_rates} 
if [ -f r4s.res ]; then
	rm r4s.res
fi 

if [ -f r4sOrig.res ]; then
	rm r4sOrig.res
fi 

if [ -f TheTree.txt ]; then
	rm TheTree.txt
fi