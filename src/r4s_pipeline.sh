#!/bin/bash
bias=$1
br_len=$2
taxa_num=$3
rep_num=$4

aln=rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.fasta
tree=n${taxa_num}_bl${br_len}.tre
r4s_norm_rates=r4s_norm_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.txt
r4s_orig_rates=r4s_orig_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.txt

##convert an alignment from nuc to aa 
python ../../../home1/02159/ds29583/r4s_benchmark/src/translate_aln.py ./r4s_benchmark_data/aln/mech_codon/nuc/${aln} ./r4s_benchmark_data/aln/mech_codon/aa/${aln}

##run rate4site 
../../../home1/02159/ds29583/rate4site.3.2.source/sourceMar09/rate4site -s ./r4s_benchmark_data/aln/mech_codon/aa/${aln} -t ./r4s_benchmark_data/trees/${tree} -o ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/r4s_rates/raw_rates/${r4s_norm_rates} -y ../../../home1/02159/ds29583/r4s_benchmark/mech_codon/r4s_rates/raw_rates/${r4s_orig_rates} 
if [ -f r4s.res ]; then
	rm r4s.res 
fi 

if [ -f r4sOrig.res ]; then
	rm r4sOrig.res  
fi 

if [ -f TheTree.txt ]; then
	rm TheTree.txt 
fi