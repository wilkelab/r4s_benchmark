#!/bin/bash
num_sim=30
site_dupl=100000
total_sites=12

if [ ! -d "$HOME/substitution_matrices_in_pheno_models/aln" ]; then
	mkdir $HOME/substitution_matrices_in_pheno_models/aln
fi

if [ ! -d "$HOME/substitution_matrices_in_pheno_models/aln/ten_sites" ]; then
	mkdir $HOME/substitution_matrices_in_pheno_models/aln/ten_sites
fi

if [ -f "$SCRATCH/run_sim_aln_ten_sites.sh" ]; then
	rm $SCRATCH/run_sim_aln_ten_sites.sh
fi


for br_len in `seq 0.02 0.02 0.50` 
do	
	for i in $(seq 1 $num_sim) 
	do
		tree_file=n2_bl${br_len}.tre
		aln_file=n2_bl${br_len}_${i}.fa
		echo "python $HOME/simulate_aln.py $HOME/substitution_matrices_in_pheno_models/132L_A_foldx_ddG.txt $HOME/substitution_matrices_in_pheno_models/trees/${tree_file} $HOME/substitution_matrices_in_pheno_models/aln/ten_sites/${aln_file} ${site_dupl} ${total_sites}" >> $SCRATCH/run_sim_aln_ten_sites.sh	
	done
done

chmod +x $SCRATCH/run_sim_aln_ten_sites.sh