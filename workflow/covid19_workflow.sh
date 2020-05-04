#!/usr/bin/env bash
SINGIMG="/srv/rbd/covid19/thecontainer/covid19_latest.sif"
GIT_PATH="$(dirname "$(readlink -f "$0")")"
#GIT_PATH=/srv/rbd/tym/covid19/workflow

runID=$1
SCHEME=$2

THREADS=100

# This is the full workflow script.
# Make final_output to dump important stuff.
mkdir -p $runID/final_output

###############################################################################
# Run demultiplexing.
###############################################################################

#singularity exec $SINGIMG bash -B /srv/rbd:/srv/rbd -c "$GIT_PATH/demultiplex.sh $runID/fastq $runID/*.csv $runID/demultiplexed $THREADS"
if [ -d $runID/demultiplexed ]; then 
  echo "demultiplexed data exists, skipping this part."
else 
  bash $GIT_PATH/demultiplex.sh $runID/fastq $runID/*.csv $runID/demultiplexed $THREADS
fi

###############################################################################
# Run processing
###############################################################################

source activate artic-ncov2019-medaka

# Medaka cannot control mem, need to scale down.
THREADS_MEDAKA=$((($THREADS+1)/3))

# Rerun artic only if not existing.
if [ -d $runID/processing/articminion ]; then
  Flag="-a"
fi
bash $GIT_PATH/processing.sh -d $runID -s nCoV-2019/$SCHEME -o $runID/processing -t $THREADS_MEDAKA $Flag
retn_code=$?

if [ $retn_code == 1 ]; then echo "ERROR in processing.sh, exitting."; exit; fi

source deactivate

###############################################################################
# Run QC
###############################################################################

source activate nextstrain

QC=$GIT_PATH/QC.sh
RMD=$GIT_PATH/QC.rmd

bash $QC -b $runID -r $RMD -t $THREADS

source deactivate

################################################################################
# Sweep important data and put in "output"
################################################################################

cp $runID/QC/$runID.html $runID/final_output/$runID.html
cp $runID/processing/results/consensus.fasta $runID/final_output/consensus.fasta
