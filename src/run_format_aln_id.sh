#!/bin/bash

hiv_aln_files=$HOME/r4s_benchmark/natural_prot/aln/hiv1/*_clean_dna_aligned.fasta
gpcr_aln_files=$HOME/r4s_benchmark/natural_prot/aln/gpcr/*_n_aligned.fasta
nuc_aln_files=( ${hiv_aln_files[@]} ${gpcr_aln_files[@]} )

for aln in ${nuc_aln_files[*]}
do
	prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
	out_aln=$HOME/r4s_benchmark/natural_prot/aln/hiv1/${prot_name}_clean_dna_reformatted.fasta
	echo $out_aln
	python src/format_aln_id.py $aln $out_aln
done