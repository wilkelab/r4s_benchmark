#!/bin/bash
rep_num=50
#bias_arr=("bias" "nobias")
bias_arr=("nobias")
taxa_num=8
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)
model='mech_codon'
distr='gamma'

if [ -f ./src/run_translate.sh ]; then
	rm ./src/run_translate.sh  
fi

for bias in ${bias_arr[*]}
do	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 4 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
				if [ $model = 'mut_sel' ]; then
					aln=rep${j}_n${i}_bl${br_len}_unequalpi_${bias}.fasta
				else
					aln=rep${j}_n${i}_bl${br_len}_${bias}.fasta
				fi
				
				if [ $distr = 'gamma' ]; then
					echo "python ./src/translate_aln_codon_to_aa.py $SCRATCH/r4s_benchmark_data/aln/mech_codon/nuc/gamma_distr/${aln} $SCRATCH/r4s_benchmark_data/aln/mech_codon/aa/gamma_distr/${aln}" >> ./src/run_translate.sh 
				else
					echo "python ./src/translate_aln_codon_to_aa.py $SCRATCH/r4s_benchmark_data/aln/${model}/nuc/${aln} $SCRATCH/r4s_benchmark_data/aln/${model}/aa/${aln}" >> ./src/run_translate.sh 
				fi
			done
		done
	done
done

chmod +x ./src/run_translate.sh