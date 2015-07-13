#!/bin/bash
sim_num=30
sim_model_arr=(dN dN_dS ms_dS ms_no_dS)
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1)

for model in ${sim_model_arr[*]}
do
	rm ./src/run_${model}_sim_pipeline.sh
	
	for num_taxa in ${taxa_num_arr[*]}
	do	
		for br_len in ${br_len_arr[*]} 
		do	 
			for i in $(seq 1 $sim_num) 
			do
    			echo "./src/sim_pipeline.sh $model $num_taxa $br_len $i" >> ./src/run_${model}_sim_pipeline.sh
			done
		done
	done
	
	chmod +x ./src/run_${model}_sim_pipeline.sh
done
