#!/bin/bash
##Setting up the parameters to make file names for simulation output
aln_file_arr=($SCRATCH/r4s_benchmark_data/aln/mech_codon/aa/gamma_distr/*)
distr='gamma'

if [ -f ./src/run_find_unchanged_sites.sh ]; then
	rm ./src/run_find_unchanged_sites.sh
fi

if [ $distr = "gamma" ]; then
	for aln in ${aln_file_arr[*]} 
	do	 
		temp=`echo $aln | grep -oP 'rep.*\.fasta$'`
		out_file=`echo $temp | sed s/.fasta/_unchanged_sites.txt/`
		model='mech_codon'
		
		echo python ./src/find_unchanged_sites.py ${aln} ./${model}/filtered_sites/gamma_distr/${out_file}  >> ./src/run_find_unchanged_sites.sh
	done
else
	for aln in ${aln_file_arr[*]} 
	do	 
		temp=`echo $aln | grep -oP 'rep.*\.fasta$'`
		out_file=`echo $temp | sed s/.fasta/_unchanged_sites.txt/`
		model='mech_codon'
		
		echo python ./src/find_unchanged_sites.py ${aln} ./${model}/filtered_sites/${out_file}  >> ./src/run_find_unchanged_sites.sh
	done
fi

natural_aln=false
if [ $natural_aln ]; then 
	hiv_aln_files=./natural_prot/aln/aligned_seqs/*_clean_protein.fasta
	gpcr_aln_files=./natural_prot/aln/aligned_seqs/*_p.fasta
	prot_aln_files=( ${hiv_aln_files[@]} ${gpcr_aln_files[@]} )

	for aln in ${prot_aln_files[*]}
	do	
		prot_name=`echo $aln | grep -oE 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
		out_file=${prot_name}_unchanged_sites.txt
		echo python ./src/find_unchanged_sites.py ${aln} ./natural_prot/filtered_sites/${out_file}  >> ./src/run_find_unchanged_sites.sh
	done
fi

chmod +x ./src/run_find_unchanged_sites.sh