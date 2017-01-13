hiv_aln_files=~/r4s_benchmark/natural_prot/aln/aligned_seqs/*_clean_protein.fasta
gpcr_aln_files=~/r4s_benchmark/natural_prot/aln/aligned_seqs/*_p.fasta
prot_aln_files=( ${hiv_aln_files[@]} ${gpcr_aln_files[@]} )

rm ./src/run_raxml.sh

for aln in ${prot_aln_files[*]}
do
	prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
	out_tree=${prot_name}.tre
	out_dir=~/r4s_benchmark/natural_prot/trees
	if [ -f ${out_dir}/RAxML_bestTree.${out_tree} ]; then
		echo ${out_dir}/RAxML_bestTree.${out_tree}
		continue 
	else
		echo $HOME/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 48 -s $aln -w $out_dir -n $out_tree -m PROTCATLG -p 1 >> ./src/run_raxml.sh
	fi
done

chmod +x ./src/run_raxml.sh
