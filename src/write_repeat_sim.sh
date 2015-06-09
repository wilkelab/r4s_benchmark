#!/bin/bash
sim_model_arr=(dN) ## dN_dS ms_dS ms_no_dS)
num_sim=30

if [ -f "./src/repeat_sim.sh" ]; then
	rm ./src/repeat_sim.sh 
fi

for m in ${sim_model_arr[*]}
do
	for i in $(seq 1 $num_sim)
	do	
		echo "./src/run_sim_pipeline.sh $m $i" >> ./src/repeat_sim.sh 
	done
done

chmod +x ./src/repeat_sim.sh