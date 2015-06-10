#!/bin/bash
sim_model_arr=(dN dN_dS) ## ms_dS ms_no_dS)
num_sim=30

if [ -f "./src/repeat_sim.sh" ]; then
	rm ./src/repeat_sim.sh 
fi

for model in ${sim_model_arr[*]}
do
	for i in $(seq 1 $num_sim)
	do	
		echo "./src/run_sim_pipeline.sh $model $i" >> ./src/repeat_sim.sh 
	done
	
	echo "Rscript ./src/plot_r4s_rates_v_sim_rates.r $model" >> ./src/repeat_sim.sh 
	
done

chmod +x ./src/repeat_sim.sh