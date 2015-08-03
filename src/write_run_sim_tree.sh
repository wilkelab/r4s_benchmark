#!/bin/bash
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1 0.33 1.0 3.3)

if [ -f "./src/run_sim_tree.sh" ]; then
	rm ./src/run_sim_tree.sh 
fi

if [ ! -d "../r4s_benchmark_data/trees" ]; then
	mkdir "../r4s_benchmark_data/trees"
fi

for num_taxa in ${taxa_num_arr[*]}
do	
	for br_len in ${br_len_arr[*]}  
	do
		tree=t${num_taxa}_b${br_len}.tre ##tree file nam
		if [ ! -d "../r4s_benchmark_data/trees/${tree}" ]; then
			echo "Rscript ./src/generate_balanced_tree_ape.R $num_taxa $br_len ../r4s_benchmark_data/trees/${tree}" >> ./src/run_sim_tree.sh
		fi
	done
done


chmod +x ./src/run_sim_tree.sh