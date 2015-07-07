#!/bin/bash
sim_num=30
sim_model_arr=(dN dN_dS ms_dS ms_no_dS)
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1)

if [ -f "./src/run_sim_tree.sh" ]; then
	rm ./src/run_sim_tree.sh 
fi

for model in ${sim_model_arr[*]}
do
if [ ! -d "${model}/trees" ]; then
	mkdir "$model"
	mkdir "${model}/trees"
fi
	for i in $(seq 1 $sim_num)
	do	
		for num_taxa in ${taxa_num_arr[*]}
		do	
			for br_len in ${br_len_arr[*]}  
			do
    			echo "Rscript src/generate_balanced_tree_ape.R $num_taxa $br_len trees/${tree}" >> ./src/run_sim_tree.sh
			done
		done
	done
done

chmod +x ./src/run_sim_tree.sh