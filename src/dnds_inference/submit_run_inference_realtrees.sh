#! /bin/bash

# SJS. Generate submission lines for dN/dS inference on Stampede. The stdout is sent to the file inference_commands


SBATCH_RAW=raw_launcher_inference.slurm

#rm launcher_rep*
#rm inference_commands*

# create output directory as needed
REPODIR=$HOME/dnmu
OUTDIR=${REPODIR}/results_realtrees
mkdir -p $OUTDIR

LAUNCHFILE=launcher1.slurm
PARAMFILE=inference_commands1
touch $LAUNCHFILE
touch $PARAMFILE
    
for TREE in h3
do
for REP in {36..50}
do    
    for TYPE in gtr bias_gtr
    do  
        for METHOD in FUBAR1 FUBAR2 SLAC_GTR                     
        do   
            DATA=rep${REP}_${TREE}_${TYPE}
            TREEFILE=$REPODIR/data/trees/${TREE}.tre 
            ALN=$REPODIR/data/alignments_realtrees/${DATA}.fasta
            OUTFILE1=$OUTDIR/${DATA}_${METHOD}.txt        # dnds inference
            OUTFILE2=$OUTDIR/${DATA}_${METHOD}_nucfit.txt # nucleotide fit

            echo sh run_inference.sh $REPODIR $DATA $ALN $TREEFILE $METHOD $OUTFILE1 $OUTFILE2 >> $PARAMFILE
        done       
    done
done
done
sed "s/PLACEHOLDER/$PARAMFILE/" $SBATCH_RAW > $LAUNCHFILE
sbatch $LAUNCHFILE

	
