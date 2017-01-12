#!/bin/bash
aln_files=./natural_prot/aln/raw_aln/*_p.fasta

if [ -f ./src/align_natural_prot.sh ]; then
	rm ./src/align_natural_prot.sh  
fi

for aln in ${aln_files[*]}
do
	out_dir=aligned_seqs
	out_file="${aln/raw_aln/$out_dir}"
	echo mafft $aln \> $out_file >> ./src/align_natural_prot.sh
done

chmod +x ./src/align_natural_prot.sh