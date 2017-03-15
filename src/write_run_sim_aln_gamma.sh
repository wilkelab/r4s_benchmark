#!/bin/bash
rep_num=50
bias_arr=("nobias")
taxa_num=9
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)
model="mech_codon"
distr="gamma"

for i in $(seq 4 $taxa_num) 
do
	param_file=run_sim_aln_n${i}.sh
				
	if [ -f $SCRATCH/$param_file ]; then
		rm $SCRATCH/$param_file
	fi

	for bias in ${bias_arr[*]}
	do	 
		for br_len in ${br_len_arr[*]} 
		do	 
			for j in $(seq 1 $rep_num) 
			do
				aln=rep${j}_n${i}_bl${br_len}_${bias}.fasta
				tree=n${i}_bl${br_len}.tre
				rates_info_file=sim_rates_info_rep${j}_n${i}_bl${br_len}_${bias}.txt
				rates_files=sim_rates_rep${j}_n${i}_bl${br_len}_${bias}.txt
				
				echo python $HOME/r4s_benchmark/src/simulate_dNdS_aln.py $bias $SCRATCH/r4s_benchmark_data/trees/$tree $SCRATCH/r4s_benchmark_data/aln/mech_codon/nuc/gamma_distr/$aln $HOME/r4s_benchmark/mech_codon/assigned_rates/raw_rates/gamma_distr/$rates_files $HOME/r4s_benchmark/mech_codon/assigned_rates/raw_rates/gamma_distr/$rates_info_file $distr >> $SCRATCH/$param_file
			done
			
		done
		
	done
	chmod +x $SCRATCH/$param_file
done

