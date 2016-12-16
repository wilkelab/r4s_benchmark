#!/bin/bash
hiv_aln_files=$HOME/r4s_benchmark/natural_prot/aln/hiv1/*_clean_protein.fasta
gpcr_aln_files=$HOME/r4s_benchmark/natural_prot/aln/gpcr/*_p_aligned.fasta
prot_aln_files=( ${hiv_aln_files[@]} ${gpcr_aln_files[@]} )

if [ -f "$SCRATCH/run_r4s_natural_prot.sh" ]; then
	rm $SCRATCH/run_r4s_natural_prot.sh
fi

for aln in ${prot_aln_files[*]}
do	
	prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
	r4s_norm_rates=$HOME/r4s_benchmark/natural_prot/r4s_rates/raw_rates/${prot_name}_norm_rates.txt
	r4s_orig_rates=$HOME/r4s_benchmark/natural_prot/r4s_rates/raw_rates/${prot_name}_orig_rates.txt
	tree_file=$HOME/r4s_benchmark/natural_prot/trees/RAxML_bestTree.${prot_name}_p.tre
	echo $HOME/rate4site.3.2.source/sourceMar09/rate4site -s $aln -t $tree_file -o $r4s_norm_rates -y $r4s_orig_rates >> $SCRATCH/run_r4s_natural_prot.sh
done

chmod +x $SCRATCH/run_r4s_natural_prot.sh