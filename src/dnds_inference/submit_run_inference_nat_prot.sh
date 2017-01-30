#! /bin/bash

# SJS. Generate submission lines for dN/dS inference on Stampede. The stdout is sent to the file inference_commands


TYPE="bias"
SBATCH_RAW=raw_launcher_inference.slurm

rm launcher_*
rm inference_commands*

# create output directory as needed
REPODIR=$HOME/r4s_benchmark/natural_prot
OUTDIR=${REPODIR}/inferred_dNdS
mkdir -p $OUTDIR

nuc_aln_files=$HOME/r4s_benchmark/natural_prot/aln/reformatted_aln/*.fasta

for aln in ${nuc_aln_files[*]}
do
    # touch launcher file and paramfile to go along with it for this rep 
    prot_name=`echo $aln | grep -oP 'hiv1_[a-zA-Z0-9]+|ENST[0-9]+'`
    
    LAUNCHFILE=launcher_${prot_name}.slurm
    PARAMFILE=inference_commands_${prot_name}
	touch $LAUNCHFILE
    touch $PARAMFILE
    METHOD=FEL1
    
    DATA=${prot_name}
    TREE=$HOME/r4s_benchmark/natural_prot/trees/reformatted_trees/${DATA}.tre
    echo $TREE
    echo $aln
    
	ALN=$aln
	OUTFILE1=$OUTDIR/${DATA}_${METHOD}.txt        # dnds inference
	OUTFILE2=$OUTDIR/${DATA}_${METHOD}_nucfit.txt # nucleotide fit
    echo sh run_inference.sh $REPODIR $DATA $ALN $TREE $METHOD $OUTFILE1 $OUTFILE2 >> $PARAMFILE
    
    # Create launcher file for this paramfile and submit. 
    sed "s/PLACEHOLDER/$PARAMFILE/" $SBATCH_RAW > $LAUNCHFILE
    sbatch $LAUNCHFILE
done

