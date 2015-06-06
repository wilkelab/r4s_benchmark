#!/bin/bash
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1)

rm run_sim.sh
for num in ${taxa_num_arr[*]}
do	
	for len in ${br_len_arr[*]}  
	do
    	echo "./sim_pipeline.sh $num $len" >> run_sim.sh
	done
done
chmod +x run_sim.sh