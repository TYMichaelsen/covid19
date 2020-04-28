#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d dir -r file] 
-- COVID-19 QC script v. $VERSION:  

Arguments:
    -h  Show this help text.
    -d  Directory with batches from 'processing.sh' 
    -b  What batch to do QC for.
    -r  R script to run to generate QC report.
    -o  Output folder.

Output:
    To come.
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hd:b:r:o:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) PROCDIR=$OPTARG;;
    b) BATCH=$OPTARG;;
    r) RSCRIPT=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${PROCDIR+x} ]; then echo "-d $MISSING"; exit 1; fi;
if [ -z ${BATCH+x} ]; then echo "-b $MISSING"; exit 1; fi;
if [ -z ${RSCRIPT+x} ]; then echo "-r $MISSING"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then echo "-o $MISSING"; exit 1; fi;

### Code.----------------------------------------------------------------------
#BATCHDIR=/srv/rbd/covid19/artic-analysis
#RSCRIPT=/srv/rbd/tym/covid19/workflow/QC.rmd
#OUTDIR=/srv/rbd/test_workflow/QC/CJ024

# Setup dirs.
mkdir -p QC
mkdir -p QC/aligntree
mkdir -p QC/$OUTDIR
mkdir -p QC/$OUTDIR/tmpdir

# Make alignment and tree of all sequences.
source activate nextstrain

bash /srv/rbd/tym/covid19/workflow/nextstrain.sh -s <(cat $PROCDIR/CJ*/results/consensus.fasta) -o QC/aligntree -t 100

source deactivate 

# Render report. 
# run QC.rmd with input args: $BATCH, $PROCDIR, $OUTDIR

exit 1

######################
# MISC STUFF


# Find lab metadata through log-file.
BATCHES=( $BATCHDIR/CJ* )

for i in "${BATCHES[@]}"; do
  METADIR=$(grep "\-d" $i/log.out | sed 's/-d: //')
  if [ -f $METADIR/*sequencing.csv ]; then 
    echo "Found
done > $OUTDIR/tmpdir/lab_meta





for i in "${BATCHES[@]}"; do grep "\-d" $i/log.out | sed 's/-d: //' ; done

grep "-d"   
cp $POOLDIR/*sequencing.csv $OUTDIR/