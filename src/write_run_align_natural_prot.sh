#!/bin/bash
prot_arr=("ENST00000000412" "ENST00000009530" "ENST00000011653" "ENST00000014914" "ENST00000023897" "ENST00000053243")

if [ -f ./src/align_natural_prot.sh ]; then
	rm ./src/align_natural_prot.sh  
fi

for prot in ${prot_arr[*]}
do
	aligned_p_file=natural_prot/aln/gpcr/${prot}_p_aligned.fasta
	aligned_n_file=natural_prot/aln/gpcr/${prot}_n_aligned.fasta
	
	if [ -f $aligned_p_file ]; then
		continue 
	else
		echo mafft natural_prot/aln/gpcr/${prot}_p.fasta \> $aligned_p_file >> ./src/align_natural_prot.sh
	fi
	
	if [ -f $aligned_n_file ]; then
		continue 
	else
		echo mafft natural_prot/aln/gpcr/${prot}_n.fasta \> $aligned_n_file >> ./src/align_natural_prot.sh
	fi
	
done

chmod +x ./src/align_natural_prot.sh