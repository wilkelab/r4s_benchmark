#!/bin/bash

nuc_aln_files=$HOME/r4s_benchmark/natural_prot/aln/back_translated_aln/*.fasta

for aln in ${nuc_aln_files[*]}
do
	prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
	out_file=$HOME/r4s_benchmark/natural_prot/aln/reformatted_aln/${prot_name}.fasta
	python src/format_aln_id.py $aln $out_file
done