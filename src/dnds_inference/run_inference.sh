#! /bin/bash

# SJS. This script runs an omega inference on Stampede.

REPODIR=$1
DATA=$2
ALN=$3
TREE=$4
METHOD=$5
OUTFILEDNDS=$6
OUTFILENUC=$7


WDIR=$SCRATCH/$DATA-$METHOD
mkdir $WDIR

BATCHFILE=auto${METHOD}.bf

cp -r batchfiles $WDIR
cp $BATCHFILE $WDIR

cd $WDIR
cp -r batchfiles/FUBAR* .
cp $ALN temp.fasta
cp $TREE temp.tre
cp $HOME/hyphy/HYPHYMP .

TOPDIR=`pwd`
sed -i "s@placeholder@$TOPDIR@g" $BATCHFILE

$HOME/hyphy/HYPHYMP CPU=48 $BATCHFILE

# Parse the output
if [[ $METHOD == "SLAC" || $METHOD == "SLAC_GTR" ]]; then
    #cp $HYPHY_CODE/parse_slac.py . # this script adds p-values into slac output file
    #python parse_slac.py $OUTFILEDNDS
    sed -i "s/ //g" slac.txt
    cp slac.txt $OUTFILEDNDS      # remove spaces, so that R can read it in properly
    cp nuc.fit $OUTFILENUC
fi

if [[ $METHOD == "FUBAR1" || $METHOD == "FUBAR2" ]]; then
    cp temp.tre.fubar.csv $OUTFILEDNDS
    cp temp.tre.gtr_fit $OUTFILENUC
fi

if [[ $METHOD == "FEL1" || $METHOD == "FEL2" || $METHOD == "FEL1_GTR" || $METHOD == "FEL2_GTR" ]]; then
    cp fel.txt $OUTFILEDNDS
    cp nuc.fit $OUTFILENUC
fi

cd
rm -r $WDIR
