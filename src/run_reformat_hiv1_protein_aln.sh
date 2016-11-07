#!/bin/bash

hiv_aln_files=$HOME/r4s_benchmark/natural_prot/aln/hiv1/*_clean_protein.fasta
 
for aln in ${hiv_aln_files[*]}
do
	prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
	out_aln=$HOME/r4s_benchmark/natural_prot/aln/hiv1/${prot_name}_clean_protein_reformatted.fasta
	echo $out_aln
	python src/format_hiv1_protein_aln.py $aln $out_aln
done