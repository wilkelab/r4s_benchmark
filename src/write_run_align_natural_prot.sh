#!/bin/bash
gpcr_aln_files=./natural_prot/aln/r4s_benchmark/natural_prot/aln/raw_aln/*_p.fasta

if [ -f ./src/align_natural_prot.sh ]; then
	rm ./src/align_natural_prot.sh  
fi

for aln in ${aln_files[*]}
do
	out_dir=aligned_seqs
	out_file="${aln/raw_aln/$out_dir}"
	echo $out_file
	#echo mafft $aln \> $aligned_n_file >> ./src/align_natural_prot.sh
done

chmod +x ./src/align_natural_prot.sh