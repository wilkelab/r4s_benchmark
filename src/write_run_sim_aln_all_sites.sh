#!/bin/bash
num_sim=30
site_dupl=1
total_sites=130
taxa_num_arr=(64 128 256 512)
br_len_arr=(0.00005 0.0005 0.005 0.05)

#if [ ! -d "$HOME/substitution_matrices_in_pheno_models/aln" ]; then
#	mkdir $HOME/substitution_matrices_in_pheno_models/aln
#fi

#if [ ! -d "$HOME/substitution_matrices_in_pheno_models/aln/all_sites" ]; then
#	mkdir $HOME/substitution_matrices_in_pheno_models/aln/all_sites
#fi

if [ -f "$SCRATCH/run_sim_aln_all_sites.sh" ]; then
	rm $SCRATCH/run_sim_aln_all_sites.sh
fi


for br_len in ${br_len_arr[*]} 
do	 
	for num in ${taxa_num_arr[*]}  
	do
		for i in $(seq 1 $num_sim) 
		do
			tree_file=n${num}_bl${br_len}.tre
			aln_file=n${num}_bl${br_len}_${i}.fa
			echo "python $HOME/substitution_matrices_in_pheno_models/src/simulate_aln.py $HOME/substitution_matrices_in_pheno_models/132L_A_foldx_ddG.txt $HOME/substitution_matrices_in_pheno_models/trees/${tree_file} $SCRATCH/substitution_matrices_in_pheno_models_data/aln/all_sites/${aln_file} ${site_dupl} ${total_sites}" >> $SCRATCH/run_sim_aln_all_sites.sh	
		done
	done
done

chmod +x $SCRATCH/run_sim_aln_all_sites.sh