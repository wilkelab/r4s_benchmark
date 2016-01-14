#!/bin/bash
rep_num=50
bias_arr=("bias" "nobias")
taxa_num=11
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

for bias in ${bias_arr[*]}
do
	if [ ! -d "$SCRATCH/r4s_benchmark_data/aln/mech_codon/nuc" ]; then
		mkdir "$SCRATCH/r4s_benchmark_data/aln/mech_codon/nuc"
	fi
		
	if [ ! -d "$SCRATCH/r4s_benchmark_data/aln/mech_codon/ancestral_aln/" ]; then
		mkdir "$SCRATCH/r4s_benchmark_data/aln/mech_codon/ancestral_aln/"
	fi	
	
	if [ ! -d "./mech_codon" ]; then
		mkdir "./mech_codon"
		mkdir "./mech_codon/true_rates"
	fi
	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 7 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
				aln=rep${j}_n${i}_bl${br_len}_${bias}.fasta
				tree=n${i}_bl${br_len}.tre
				true_rates=sim_rates_rep${j}_n${i}_bl${br_len}_${bias}.txt
				echo "python ../../../home1/02159/ds29583/r4s_benchmark/src/simulate_mech_codon_aln.py $bias ./r4s_benchmark_data/trees/${tree} ./r4s_benchmark_data/aln/mech_codon/nuc/$aln ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/true_rates/${true_rates}" >> ./src/run_sim_aln.sh
			done
		done
	done
done

chmod +x ./src/run_sim_aln.sh