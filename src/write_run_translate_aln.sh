#!/bin/bash
rep_num=50
bias_arr=("bias" "nobias")
taxa_num=11
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)
model='mut_sel'

if [ -f ./src/run_translate.sh ]; then
	rm ./src/run_translate.sh  
fi

for bias in ${bias_arr[*]}
do	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 7 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
				if [ $model = 'mut_sel' ]; then
					aln=rep${j}_n${i}_bl${br_len}_unequalpi_${bias}.fasta
				else
					aln=rep${j}_n${i}_bl${br_len}_${bias}.fasta
				fi
				
				echo "python ../../../home1/02159/ds29583/r4s_benchmark/src/translate_aln.py ./r4s_benchmark_data/aln/${model}/nuc/${aln} ./r4s_benchmark_data/aln/${model}/aa/${aln}" >> ./src/run_translate.sh 
			done
		done
	done
done

chmod +x ./src/run_translate.sh