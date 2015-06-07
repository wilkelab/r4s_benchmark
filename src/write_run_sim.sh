#!/bin/bash
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1)

if [ -f "src/run_sim.sh" ]; then
	rm src/run_sim.sh
fi

for num in ${taxa_num_arr[*]}
do	
	for len in ${br_len_arr[*]}  
	do
    	echo "./src/sim_pipeline.sh $num $len" >> src/run_sim.sh
	done
done
chmod +x ./src/run_sim.sh