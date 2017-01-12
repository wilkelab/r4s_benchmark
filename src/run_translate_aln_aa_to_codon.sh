#! /bin/bash

hiv_aln_files=natural_prot/aln/aligned_seqs/*_clean_protein.fasta
gpcr_aln_files=natural_prot/aln/aligned_seqs/*_p.fasta

for aln in ${gpcr_aln_files[*]}
do
	prot_name=`echo $aln | grep -o -E 'ENST[0-9]+'`
	codon_file=${prot_name}_n.fasta
	python src/translate_aln_aa_to_codon.py $aln natural_prot/aln/raw_aln/$codon_file natural_prot/aln/back_translated_aln/$codon_file
done

for aln in ${hiv_aln_files[*]}
do
	prot_name=`echo $aln | grep -o -E 'hiv1_[a-zA-Z0-9]+'`
	codon_file=${prot_name}_clean_dna.fasta
	python src/translate_aln_aa_to_codon.py $aln natural_prot/aln/raw_aln/$codon_file natural_prot/aln/back_translated_aln/$codon_file
done