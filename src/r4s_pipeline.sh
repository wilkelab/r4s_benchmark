#!/bin/bash
bias=$1
br_len=$2
taxa_num=$3
rep_num=$4
model=$5
distr=$6
gamma_num=$7

if [ $model = 'mut_sel' ]; then
	aln=rep${rep_num}_n${taxa_num}_bl${br_len}_unequalpi_${bias}.fasta
else
	aln=rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.fasta
fi

tree=n${taxa_num}_bl${br_len}.tre
tree_dir=./r4s_benchmark_data/trees

if [ $distr = 'gamma' ]; then
	aln=rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}_gamma${gamma_num}.fasta
	aln_dir=./r4s_benchmark_data/aln/${model}/aa/gamma_distr
	rates_dir=../../../home1/02159/ds29583/r4s_benchmark/${model}/r4s_rates/raw_rates/gamma_distr
	r4s_norm_rates=r4s_norm_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}_gamma${gamma_num}.txt
	r4s_orig_rates=r4s_orig_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}_gamma${gamma_num}.txt
else
	aln_dir=./r4s_benchmark_data/aln/${model}/aa
	rates_dir=../../../home1/02159/ds29583/r4s_benchmark/${model}/r4s_rates/raw_rates
	r4s_norm_rates=r4s_norm_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.txt
	r4s_orig_rates=r4s_orig_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.txt
fi

##run rate4site 
../../../home1/02159/ds29583/rate4site.3.2.source/sourceMar09/rate4site -s $aln_dir/$aln -t $tree_dir/$tree -o $rates_dir/$r4s_norm_rates -y $rates_dir/$r4s_orig_rates
if [ -f r4s.res ]; then
	rm r4s.res 
fi 

if [ -f r4sOrig.res ]; then
	rm r4sOrig.res  
fi 

if [ -f TheTree.txt ]; then
	rm TheTree.txt 
fi