#!/bin/bash
##Setting up the parameters to make file names for simulation output
rep_num=50
bias_arr=("bias" "nobias")
taxa_num=11
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)

##Creating directory structure for the ../r4s_benchmark_data and mech_codon
if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

if [ ! -d "$SCRATCH/r4s_benchmark_data/aln/mech_codon/nuc" ]; then
	mkdir $SCRATCH/r4s_benchmark_data/aln/mech_codon/nuc
fi

if [ ! -d "$SCRATCH/r4s_benchmark_data/aln/mech_codon/ancestral" ]; then
	mkdir $SCRATCH/r4s_benchmark_data/aln/mech_codon/ancestral
fi	

if [ ! -d "mech_codon/assigned_rates/" ]; then
	mkdir mech_codon/assigned_rates/
fi

if [ ! -d "mech_codon/assigned_rates/raw_rates" ]; then
	mkdir mech_codon/assigned_rates/raw_rates
fi

if [ ! -d "mech_codon/assigned_rates/processed_rates" ]; then
	mkdir mech_codon/assigned_rates/processed_rates
fi

for bias in ${bias_arr[*]}
do	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 7 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
				aln=rep${j}_n${i}_bl${br_len}_${bias}.fasta
				tree=n${i}_bl${br_len}.tre
				sim_rates=sim_rates_rep${j}_n${i}_bl${br_len}_${bias}.txt
				sim_rates_info=sim_rates_info_rep${j}_n${i}_bl${br_len}_${bias}.txt
				true_rates=true_rates_rep${j}_n${i}_bl${br_len}_${bias}.txt
				echo "python ../../../home1/02159/ds29583/r4s_benchmark/src/simulate_mech_codon_aln.py $bias ./r4s_benchmark_data/trees/${tree} ./r4s_benchmark_data/aln/mech_codon/nuc/${aln} ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/assigned_rates/raw_rates/${sim_rates} ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/assigned_rates/raw_rates/${sim_rates_info}"  >> ./src/run_sim_aln.sh
			done
		done
	done
done

chmod +x ./src/run_sim_aln.sh