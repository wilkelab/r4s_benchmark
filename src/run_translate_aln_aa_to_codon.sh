#! /bin/bash

hiv_aln_files=natural_prot/aln/hiv1/*_clean_protein.fasta
gpcr_aln_files=natural_prot/aln/gpcr/*_p_aligned.fasta

for aln in ${gpcr_aln_files[*]}
do
	prot_name=`echo $aln | grep -o -E 'ENST[0-9]+'`
	unaln_codon_file=${prot_name}_n.fasta
	out_aln=${prot_name}_n_aligned.fasta
	python src/translate_aln_aa_to_codon.py $aln natural_prot/aln/gpcr/$unaln_codon_file natural_prot/aln/gpcr/$out_aln
done

for aln in ${hiv_aln_files[*]}
do
	prot_name=`echo $aln | grep -o -E 'hiv1_[a-zA-Z0-9]+'`
	unaln_codon_file=${prot_name}_clean_dna.fasta
	out_aln=${prot_name}_clean_dna_aligned.fasta
	python src/translate_aln_aa_to_codon.py $aln natural_prot/aln/hiv1/$unaln_codon_file natural_prot/aln/hiv1/$out_aln
done