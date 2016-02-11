#!/bin/bash
##Setting up parameters to call files
rep_num=50
bias_arr=("bias" "nobias")
taxa_num=11
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)

if [ -f "./src/run_merge.sh" ]; then
	rm ./src/run_merge.sh 
fi

##Create the necessary directory structure in the ./mech_codon for the run_merge.sh output.
if [ ! -d "mech_codon/sim_rates/assigned_rates/processed_rates" ]; then
		mkdir mech_codon/sim_site_rates/assigned_rates/processed_rates
fi
	
##Write run_merge.sh file to merge pyvolve output files sim_rates and sim_rates_info into one file that contains information on the assigned site-specific dN/dS for each site.
for bias in ${bias_arr[*]}
do	 
	for br_len in ${br_len_arr[*]} 
	do	 
		for i in $(seq 7 $taxa_num) 
		do
			for j in $(seq 1 $rep_num) 
			do
				##generate a file name using the parameters mentioned above
				sim_rates=sim_rates_rep${j}_n${i}_bl${br_len}_${bias}.txt
				sim_rates_info=sim_rates_info_rep${j}_n${i}_bl${br_len}_${bias}.txt
    			echo "Rscript ../../../home1/02159/ds29583/r4s_benchmark/src/merge_site_rates.r ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/sim_rates/assigned_rates/raw_rates/${sim_rates} ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/sim_rates/assigned_rates/raw_rates/${sim_rates_info} ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/sim_rates/assigned_rates/processed_rates/${sim_rates}" >> ./src/run_merge.sh
			done
		done
	done
done
chmod +x ./src/run_merge.sh