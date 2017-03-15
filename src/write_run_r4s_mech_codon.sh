#!/bin/bash
rep_num=50
#bias_arr=("bias" "nobias")
bias_arr=('nobias')
#taxa_num=11
taxa_num=8
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)
model='mech_codon'

if [ ! -d "../../../home1/02159/ds29583/r4s_benchmark/${model}/r4s_rates" ]; then
	mkdir ../../../home1/02159/ds29583/r4s_benchmark/${model}/r4s_rates
fi

if [ ! -d "../../../home1/02159/ds29583/r4s_benchmark/${model}/r4s_rates/raw_rates" ]; then
	mkdir ../../../home1/02159/ds29583/r4s_benchmark/${model}/r4s_rates/raw_rates
fi

if [ ! -d "$SCRATCH/r4s_benchmark_data/aln/${model}/aa/" ]; then
	mkdir $SCRATCH/r4s_benchmark_data/aln/${model}/aa/
fi

for br_len in ${br_len_arr[*]} 
do	
	rm ../../../home1/02159/ds29583/run_bl${br_len}_r4s_mech_codon.sh 
done

for bias in ${bias_arr[*]}
do	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 4 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
				echo "../../../home1/02159/ds29583/r4s_benchmark/src/r4s_pipeline.sh $bias $br_len $i $j $model" >> ../../../home1/02159/ds29583/run_bl${br_len}_r4s_mech_codon.sh
			done
		done
	chmod +x ../../../home1/02159/ds29583/run_bl${br_len}_r4s_mech_codon.sh 
	done
done


