#!/bin/bash
model=$1
sim_num=$2
taxa_num_arr=(32 64 128 256)
br_len_arr=(0.001 0.0033 0.01 0.033 0.1)

for num_taxa in ${taxa_num_arr[*]}
do	
	for br_len in ${br_len_arr[*]}  
	do
    	./src/sim_pipeline.sh $model $sim_num $num_taxa $br_len
	done
done

