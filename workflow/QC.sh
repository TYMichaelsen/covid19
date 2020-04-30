#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d dir -r file] 
-- COVID-19 QC script v. $VERSION:  

Arguments:
    -h  Show this help text.
    -b  What batch to do QC for.
    -r  R script to run to generate QC report.
    -t  Number of threads.

Output:
    To come.
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hb:r:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    b) BATCH=$OPTARG;;
    r) RMD=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${BATCH+x} ]; then echo "-b $MISSING"; exit 1; fi;
if [ -z ${RMD+x} ]; then echo "-r $MISSING"; exit 1; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------
# Setup dirs.
mkdir -p QC
mkdir -p QC/$BATCH
mkdir -p QC/$BATCH/aligntree


###############################################################################
# Setup data to be used in QC.
###############################################################################

# Copy over sequences.
cp processing/$BATCH/results/consensus.fasta QC/$BATCH/aligntree/sequences.fasta
#cat export/*_export/sequences.fasta > QC/aligntree/sequences.fasta

### Alignment ###
augur align \
--sequences QC/$BATCH/aligntree/sequences.fasta \
--reference-sequence auxdata/reference/MN908947.3.gb \
--output QC/$BATCH/aligntree/aligned.fasta \
--nthreads $THREADS &> QC/$BATCH/aligntree/log.out

### Mask bases ###
mask_sites="18529 29849 29851 29853"

python3 auxdata/ncov/scripts/mask-alignment.py \
--alignment QC/$BATCH/aligntree/aligned.fasta \
--mask-from-beginning 130 \
--mask-from-end 50 \
--mask-sites $mask_sites \
--output QC/$BATCH/aligntree/masked.fasta
    
### Tree ###
augur tree \
--alignment QC/$BATCH/aligntree/masked.fasta \
--output QC/$BATCH/aligntree/tree_raw.nwk \
--nthreads $THREADS
  
###############################################################################
# Generate the QC report.
###############################################################################

# Fetch path lab metadata.
pth=$(grep "\-d" processing/$BATCH/log.out | sed 's/-d: //')
labmeta=$(find $pth/*sequencing.csv)

if [ ! -f $labmeta ]; then 
  echo "ERROR: could not find lab metadata. Searched for '$pth/*sequencing.csv', but found nothing."
  exit 1
fi

# Run .rmd script.
Rscript -e "rmarkdown::render(input='$RMD',output_file='/srv/rbd/covid19/current/QC/$BATCH/$BATCH.html',knit_root_dir='$PWD',params=list(batch='$BATCH',labmeta='$labmeta'))"


