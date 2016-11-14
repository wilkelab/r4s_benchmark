#! /bin/bash

nuc_aln_files=$HOME/r4s_benchmark/natural_prot/aln/*/*_clean_dna_reformatted.fasta

rm ./src/run_raxml.sh

for aln in ${nuc_aln_files[*]}
do
	prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
	out_tree=${prot_name}_n.tre
	out_dir=~/r4s_benchmark/natural_prot/trees/
	echo $HOME/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 48 -s $aln -w $out_dir -n $out_tree -m GTRGAMMA -p 1 >> ./src/run_raxml.sh
done

# hiv_aln_files=$HOME/r4s_benchmark/natural_prot/aln/hiv1/*_clean_protein.fasta
# gpcr_aln_files=$HOME/r4s_benchmark/natural_prot/aln/gpcr/*_p_aligned.fasta
# prot_aln_files=( ${hiv_aln_files[@]} ${gpcr_aln_files[@]} )
# 
# for aln in ${prot_aln_files[*]}
# do
# 	prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
# 	out_tree=${prot_name}_p.tre
# 	out_dir=~/r4s_benchmark/natural_prot/trees/
# 	echo $HOME/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 48 -s $aln -w $out_dir -n $out_tree -m PROTCATLG -p 1 >> ./src/run_raxml.sh
# done
# 
# chmod +x ./src/run_raxml.sh
