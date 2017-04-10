#!/bin/bash
rep_num=50
bias="nobias"
taxa_num=8
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)
model="mech_codon"
distr="gamma"
shape_arr=(0.23 0.35 0.647 0.312 0.378 0.238)
rate_arr=(0.029 0.024 0.070 0.032 0.062 0.018)

for i in $(seq 4 $taxa_num) 
do
	param_file=run_sim_aln_n${i}.sh
				
	if [ -f $SCRATCH/$param_file ]; then
		rm $SCRATCH/$param_file
	fi

	for br_len in ${br_len_arr[*]} 
	do	 
		for j in $(seq 1 $rep_num) 
		do
			for ((k=0; k<${#shape_arr[*]}; k++));
			do
				gamma_num=$((k+1))
				shape=${shape_arr[k]}
				rate=${rate_arr[k]}
				aln=rep${j}_n${i}_bl${br_len}_${bias}_gamma${gamma_num}.fasta
				tree=n${i}_bl${br_len}.tre
				rates_info_file=sim_rates_info_rep${j}_n${i}_bl${br_len}_${bias}_gamma${gamma_num}.txt
				rates_files=sim_rates_rep${j}_n${i}_bl${br_len}_${bias}_gamma${gamma_num}.txt
			
				echo python $HOME/r4s_benchmark/src/simulate_dNdS_aln.py $bias $SCRATCH/r4s_benchmark_data/trees/$tree $SCRATCH/r4s_benchmark_data/aln/mech_codon/nuc/gamma_distr/$aln $HOME/r4s_benchmark/mech_codon/assigned_rates/raw_rates/gamma_distr/$rates_files $HOME/r4s_benchmark/mech_codon/assigned_rates/raw_rates/gamma_distr/$rates_info_file $distr $shape $rate>> $SCRATCH/$param_file
			done
		done
	done
	chmod +x $SCRATCH/$param_file
done
	


