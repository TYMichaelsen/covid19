#!/bin/bash
runID=$1

echo $runID

# This is the full workflow script.

###############################################################################
# Run basecalling  
###############################################################################

# Some script for basecalling.

###############################################################################
# Run demultiplexing
###############################################################################

# Some script for demultiplexing.

###############################################################################
# Run processing
###############################################################################

source activate artic-ncov2019-medaka

PROCESS=/srv/rbd/tym/covid19/workflow/processing.sh

bash $PROCESS -d $runID/rawdata -s nCoV-2019/V3.1 -o $runID/processing -t 100 

source deactivate

exit 1

###############################################################################
# Run QC
###############################################################################
source activate nextstrain

QC=/srv/rbd/tym/covid19/workflow/QC.sh
RMD=/srv/rbd/tym/covid19/workflow/QC.rmd

bash $QC -b $runID -r $RMD -t 100

source deactivate
