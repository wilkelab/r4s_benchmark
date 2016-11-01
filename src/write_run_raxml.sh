#! /bin/bash

hiv_aln_files=natural_prot/aln/hiv1/*_clean_dna.fasta
gpcr_aln_files=natural_prot/aln/gpcr/*_n_aligned.fasta
nuc_aln_files=( ${hiv_aln_files[@]} ${gpcr_aln_files[@]} )

rm ./src/run_raxml.sh

for aln in ${nuc_aln_files[*]}
do
	prot_name=`echo $aln | grep -o -E 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
	out_tree=${prot_name}.tre
	out_dir=~/Desktop/projects/r4s_benchmark/natural_prot/trees/
	echo ../standard-RAxML/raxmlHPC-SSE3 -s $aln -w $out_dir -n $out_tree -m GTRGAMMA -p 1 >> ./src/run_raxml.sh
done

chmod +x ./src/run_raxml.sh