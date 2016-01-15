#!/bin/bash
rep_num=50
bias_arr=("bias" "nobias")
taxa_num=11
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)

if [ ! -d "./mech_codon/r4s_rates" ]; then
	mkdir "./mech_codon/r4s_rates"
fi

if [ ! -d "./mech_codon/r4s_rates/raw_rates" ]; then
	mkdir "./mech_codon/r4s_rates/raw_rates"
fi

if [ ! -d "$SCRATCH/r4s_benchmark_data/aln/mech_codon/aa/" ]; then
	mkdir "$SCRATCH/r4s_benchmark_data/aln/mech_codon/aa/"
fi

for br_len in ${br_len_arr[*]} 
do	
	rm ./src/run_bl${br_len}_r4s.sh
done

for bias in ${bias_arr[*]}
do	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 7 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
			echo "../../../home1/02159/ds29583/r4s_benchmark/src/r4s_pipeline.sh $bias $br_len $i $j" >> ./src/run_bl${br_len}_r4s.sh
			done
		done
	chmod +x ./src/run_bl${br_len}_r4s.sh
	done
done


