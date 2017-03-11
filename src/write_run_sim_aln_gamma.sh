#!/bin/bash
rep_num=50
bias_arr=("nobias")
taxa_num=11
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)
model="mech_codon"

if [ -f "$SCRATCH/run_sim_aln.sh" ]; then
	rm $SCRATCH/run_sim_aln.sh
fi

for bias in ${bias_arr[*]}
do	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 4 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
				aln=rep${j}_n${i}_bl${br_len}_${bias}.fasta
				tree=n${i}_bl${br_len}.tre
				rates_info_file=sim_rates_info_rep${j}_n${i}_bl${br_len}_${bias}.txt
				rates_files=sim_rates_rep${j}_n${i}_bl${br_len}_${bias}.txt
				echo python src/simulate_dNdS_aln.py $bias ../r4s_benchmark_data/trees/$tree ../r4s_benchmark_data/aln/$aln_file mech_codon/assigned_rates/$rates_files mech_codon/assigned_rates/$rates_info_file > ./src/run_sim_aln.sh
			done
		done
	done
done

chmod +x ./src/run_sim_aln.sh