#!/bin/bash
sim_num=30
sim_model_arr=(dN dN_dS ms_dS ms_no_dS)
taxa_num_arr=(32 64 128 256)
br_len_arr1=(0.001 0.0033 0.01 0.033 0.1)
br_len_arr2=(0.001 0.0033 0.01 0.033 0.1 0.33 1.0 3.3)

for model in ${sim_model_arr[*]}
do
	rm ./src/run_${model}_r4s.sh
	if [ $model = "dN" -o $model = "dN_dS" ]; then 
		br_len_arr=("${br_len_arr1[*]}")
	else 
		br_len_arr=("${br_len_arr2[*]}")
	fi
	for num_taxa in ${taxa_num_arr[*]}
	do	
		for br_len in ${br_len_arr[*]} 
		do	 
			for i in $(seq 1 $sim_num) 
			do
				r4s_norm_rates=r4s_norm_rates_t${num_taxa}_b${br_len}_${i}.txt ##rate4site output file name for norm rates
				if [ ! -f "${model}/r4s_site_rates/${r4s_norm_rates}" ]; then
					echo "./src/r4s_pipeline.sh $model $num_taxa $br_len $i"
    				echo "./src/r4s_pipeline.sh $model $num_taxa $br_len $i" >> ./src/run_${model}_r4s.sh
    			fi
			done
		done
	done
	
	chmod +x ./src/run_${model}_r4s.sh
done
