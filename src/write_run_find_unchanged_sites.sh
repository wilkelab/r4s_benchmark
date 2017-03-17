#!/bin/bash
##Setting up the parameters to make file names for simulation output
rep_num=50
#bias_arr=("bias" "nobias")
bias_arr=("nobias")
#taxa_num=11
taxa_num=8
br_len_arr=(0.0025 0.01 0.04 0.16 0.64)
#model_arr=("mut_sel" "mech_codon")
model_arr=("mech_codon")
distr='gamma'

if [ -f ./src/run_find_unchanged_sites.sh ]; then
	rm ./src/run_find_unchanged_sites.sh
fi

for model in ${model_arr[*]}
do 
	for br_len in ${br_len_arr[*]} 
	do	 
		for bias in ${bias_arr[*]}
		do	 
			#for i in $(seq 7 $taxa_num)
			for i in $(seq 4 $taxa_num) 
			do
				for j in $(seq 1 $rep_num) 
				do
					if [ $model = "mut_sel" ]; then
						aln=rep${j}_n${i}_bl${br_len}_unequalpi_${bias}.fasta
					else
						aln=rep${j}_n${i}_bl${br_len}_${bias}.fasta
					fi
					
					out_file=rep${j}_n${i}_bl${br_len}_${bias}_unchanged_sites.txt
					
					if [ $distr = "gamma" ]; then
						echo "python ./src/find_unchanged_sites.py $SCRATCH/r4s_benchmark_data/aln/${model}/aa/gamma_distr/${aln} ./${model}/filtered_sites/gamma_distr/${out_file}"  >> ./src/run_find_unchanged_sites.sh
					else
						echo "python ./src/find_unchanged_sites.py ../r4s_benchmark_data/aln/${model}/aa/${aln} ./${model}/filtered_sites/${out_file}"  >> ./src/run_find_unchanged_sites.sh
					fi
				done
			done
		done
	done
done
# 
# hiv_aln_files=./natural_prot/aln/aligned_seqs/*_clean_protein.fasta
# gpcr_aln_files=./natural_prot/aln/aligned_seqs/*_p.fasta
# prot_aln_files=( ${hiv_aln_files[@]} ${gpcr_aln_files[@]} )
# 
# for aln in ${prot_aln_files[*]}
# do	
# 	prot_name=`echo $aln | grep -oE 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
# 	out_file=${prot_name}_unchanged_sites.txt
# 	echo "python ./src/find_unchanged_sites.py ${aln} ./natural_prot/filtered_sites/${out_file}"  >> ./src/run_find_unchanged_sites.sh
# done

chmod +x ./src/run_find_unchanged_sites.sh