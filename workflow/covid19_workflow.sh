#!/bin/bash
runID=$1

AAU_COVID19_PATH="$(dirname "$(readlink -f "$0")")"

#AAU_COVID19_PATH=/srv/rbd/tym/covid19/workflow

# This is the full workflow script.

THREADS=100

###############################################################################
# Run demultiplexing
###############################################################################

#bash $AAU_COVID19_PATH/demultiplex.sh $runID/fastq $runID/*.csv $runID/demultiplexed $THREADS

###############################################################################
# Run processing
###############################################################################

source activate artic-ncov2019-medaka

bash $AAU_COVID19_PATH/processing.sh -d $runID -s nCoV-2019/V3.1 -o $runID/processing -t $THREADS

source deactivate

###############################################################################
# Run QC
###############################################################################

source activate nextstrain

QC=$AAU_COVID19_PATH/QC.sh
RMD=$AAU_COVID19_PATH/QC.rmd

bash $QC -b $runID -r $RMD -t $THREADS

source deactivate

################################################################################
# Sweep important data and put in "output"
################################################################################

mkdir $runID/final_output
cp $runID/QC/$runID.html $runID/final_output/$runID.html
cp $runID/processing/results/consensus.fasta $runID/final_output/consensus.fasta
