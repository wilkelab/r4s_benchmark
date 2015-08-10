#!/bin/bash
sim_num=30
sim_model_arr=("mech_codon_dN" "mech_codon_dN_dS")
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1)

if [ -f "./src/run_merge.sh" ]; then
	rm ./src/run_merge.sh 
fi

if [ ! -d "plots" ]; then
	mkdir plots
fi

for model in ${sim_model_arr[*]}
do	
	if [ ! -d "${model}/sim_site_rates/merged_output/" ]; then
		mkdir ${model}/sim_site_rates/merged_output/
	fi

	for num_taxa in ${taxa_num_arr[*]}
	do	
		for br_len in ${br_len_arr[*]}  
		do
			for i in $(seq 1 $sim_num)
			do
				sim_rates=site_rates_t${num_taxa}_b${br_len}_${i}.txt 
				sim_rates_info=site_rates_info_t${num_taxa}_b${br_len}_${i}.txt ##pyvolve output file name
    			echo "Rscript ./src/merge_site_rates.r ${model}/sim_site_rates/assigned_rates/${sim_rates} ${model}/sim_site_rates/assigned_rates/${sim_rates_info} ${model}/sim_site_rates/merged_output/${sim_rates}" >> ./src/run_merge.sh
			done
		done
	done
done

chmod +x ./src/run_merge.sh