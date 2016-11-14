#!/bin/bash

tree_files=$HOME/r4s_benchmark/natural_prot/trees/RAxML_bestTree.*.tre

for tree in ${tree_files[*]}
do
	prot_name=`echo $tree | grep -oP '(hiv1_[a-zA-Z0-9]+|ENST[0-9]+)_[np]'`
	out_tree=$HOME/r4s_benchmark/natural_prot/trees/RAxML_bestTree.${prot_name}_reformatted.tre
	echo $out_tree
	python src/format_tree_id.py $tree $out_tree
done