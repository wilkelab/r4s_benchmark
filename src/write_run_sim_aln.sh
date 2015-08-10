#!/bin/bash
sim_model_arr=(mech_codon_dN mech_codon_dN_dS mut_sel_dN mut_sel_dN_dS)
taxa_num_arr=(32 64 128 256)
br_len_arr1=(0.001 0.0033 0.01 0.033 0.1)
br_len_arr2=(0.001 0.0033 0.01 0.033 0.1 0.33 1.0 3.3)
sim_num=30

if [ -f "./src/run_sim_aln.sh" ]; then
	rm ./src/run_sim_aln.sh 
fi

for model in ${sim_model_arr[*]}
do
	if [ ! -d "../r4s_benchmark_data/${model}" ]; then
		mkdir "../r4s_benchmark_data/${model}"
		mkdir "../r4s_benchmark_data/${model}/aln"
		mkdir "../r4s_benchmark_data/${model}/aln/nuc"
	fi
	
	if [ ! -d "$model" ]; then
		mkdir "$model"
		mkdir "${model}/sim_site_rates/"
		mkdir "${model}/sim_site_rates/assigned_rates"
		mkdir "${model}/sim_site_rates/inferred_rates"
	fi 
	
	if [ $model = "mech_codon_dN" -o $model = "mech_codon_dN_dS" ]; then 
		br_len_arr=("${br_len_arr1[*]}")
	else 
		br_len_arr=("${br_len_arr2[*]}")
	fi
	
	echo "cd $model" >> ./src/run_sim_aln.sh
	
	for num_taxa in ${taxa_num_arr[*]}
	do	
		for br_len in ${br_len_arr[*]}  
		do
			for i in $(seq 1 $sim_num)
			do	
				tree=t${num_taxa}_b${br_len}.tre ##tree file nam
				aln=seq_t${num_taxa}_b${br_len}_${i}.fasta ##multiple sequence alignment file name
				sim_rates=site_rates_t${num_taxa}_b${br_len}_${i}.txt ##pyvolve output file name
				sim_rates_info=site_rates_info_t${num_taxa}_b${br_len}_${i}.txt
				if [ ! -f "../r4s_benchmark_data/${model}/aln/nuc/${aln}" ]; then
    				echo "python ../src/simulate_aln.py $model ../../r4s_benchmark_data/trees/${tree} ../../r4s_benchmark_data/${model}/aln/nuc/$aln sim_site_rates/assigned_rates/${sim_rates} sim_site_rates/assigned_rates/${sim_rates_info}" >> ./src/run_sim_aln.sh
				fi
			done
		done
	done
	echo "cd .." >> ./src/run_sim_aln.sh
done

chmod +x ./src/run_sim_aln.sh