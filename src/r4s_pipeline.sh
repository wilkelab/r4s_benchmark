#!/bin/bash
bias=$1
br_len=$2
taxa_num=$3
rep_num=$4

aln=rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.fasta
tree=n${taxa_num}_bl${br_len}.tre
r4s_norm_rates=r4s_norm_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.txt
r4s_orig_rates=r4s_orig_rates_rep${rep_num}_n${taxa_num}_bl${br_len}_${bias}.txt

if [ ! -d "$SCRATCH/ev_rate_method_comparison/aln/aa/" ]; then
	mkdir $SCRATCH/ev_rate_method_comparison/aln/aa/
fi

##convert an alignment from nuc to aa 
python $SCRATCH/ev_rate_method_comparison/src/translate_aln.py $SCRATCH/ev_rate_method_comparison/aln/nuc/${aln} $SCRATCH/ev_rate_method_comparison/aln/aa/${aln}

if [ ! -d "$SCRATCH/ev_rate_method_comparison/r4s_site_rates/" ]; then
	mkdir $SCRATCH/ev_rate_method_comparison/r4s_site_rates/
	mkdir $SCRATCH/ev_rate_method_comparison/r4s_site_rates/${bias}
fi

##run rate4site 
rate4site.3.2.source/sourceMar09/rate4site -s $SCRATCH/ev_rate_method_comparison/aln/aa/${aln} -t $SCRATCH/ev_rate_method_comparison/trees/${tree} -o $SCRATCH/ev_rate_method_comparison/r4s_site_rates/${r4s_norm_rates} 
if [ -f r4s.res ]; then
	rm r4s.res 
fi 

if [ -f r4sOrig.res ]; then
	mv r4sOrig.res $SCRATCH/ev_rate_method_comparison/r4s_site_rates/${r4s_orig_rates} 
fi 

if [ -f TheTree.txt ]; then
	rm TheTree.txt 
fi