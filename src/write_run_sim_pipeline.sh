#!/bin/bash
sim_num=1
sim_model_arr=(dN dN_dS ms_dS ms_no_dS)
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1)

if [ -f "./src/run_sim_pipeline.sh" ]; then
	rm ./src/run_sim_pipeline.sh 
fi

for model in ${sim_model_arr[*]}
do
	for i in $(seq 1 $sim_num)
	do	
		for num_taxa in ${taxa_num_arr[*]}
		do	
			for br_len in ${br_len_arr[*]}  
			do
    			echo "./src/sim_pipeline.sh $model $i $num_taxa $br_len" >> ./src/run_sim_pipeline.sh
			done
		done
	done
	##echo "Rscript ./src/plot_r4s_rates_v_sim_rates.r $model" >> ./src/run_sim_pipeline.sh
done

chmod +x ./src/run_sim_pipeline.sh