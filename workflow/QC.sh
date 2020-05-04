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

AAU_COVID19_PATH="$(dirname "$(readlink -f "$0")")"

# Setup dirs.
mkdir -p $BATCH/QC
mkdir -p $BATCH/QC/aligntree

# Logging
LOG_NAME="$BATCH/QC/QC_log_$(date +"%Y-%m-%d-%T").txt"
echo "QC log" >> $LOG_NAME
echo "AAU COVID-19 revision - $(git -C $AAU_COVID19_PATH rev-parse --short HEAD)" >> $LOG_NAME
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

REF=$AAU_COVID19_PATH/MN908947.3.gb
CLADES=/srv/rbd/covid19/current/auxdata/root_seqs/pangolin_clades.tsv

###############################################################################
# Setup data to be used in QC.
###############################################################################

# Copy over sequences.
cp $BATCH/processing/results/consensus.fasta $BATCH/QC/aligntree/sequences.fasta
#cat export/*_export/sequences.fasta > QC/aligntree/sequences.fasta

### Alignment ###
augur align \
--sequences $BATCH/QC/aligntree/sequences.fasta \
--reference-sequence $REF \
--output $BATCH/QC/aligntree/aligned.fasta \
--nthreads $THREADS &> $BATCH/QC/aligntree/log.out

### Mask bases ###
mask_sites="18529 29849 29851 29853"

python3 $AAU_COVID19_PATH/mask-alignment.py \
--alignment $BATCH/QC/aligntree/aligned.fasta \
--mask-from-beginning 130 \
--mask-from-end 50 \
--mask-sites $mask_sites \
--output $BATCH/QC/aligntree/masked.fasta
    
### Tree ###
augur tree \
--alignment $BATCH/QC/aligntree/masked.fasta \
--output $BATCH/QC/aligntree/tree_raw.nwk \
--nthreads $THREADS
  
augur refine \
--tree $BATCH/QC/aligntree/tree_raw.nwk \
--output-tree $BATCH/QC/aligntree/tree.nwk

### ancestral tree.
augur ancestral \
  --tree $BATCH/QC/aligntree/tree.nwk \
  --alignment $BATCH/QC/aligntree/masked.fasta \
  --output-node-data $BATCH/QC/aligntree/nt_muts.json \
  --inference joint \
  --infer-ambiguous
  
### Translate NT ot AA.
augur translate \
  --tree $BATCH/QC/aligntree/tree.nwk \
  --ancestral-sequences $BATCH/QC/aligntree/nt_muts.json \
  --reference-sequence $REF \
  --output-node-data $BATCH/QC/aligntree/aa_muts.json
                       
### add clades.
augur clades \
  --tree $BATCH/QC/aligntree/tree.nwk \
  --mutations $BATCH/QC/aligntree/nt_muts.json $BATCH/QC/aligntree/aa_muts.json \
  --clades $CLADES \
  --output-node-data $BATCH/QC/aligntree/clades.json
  
###############################################################################
# Generate the QC report.
###############################################################################

# Fetch path lab metadata.
pth=$(grep "\-d" $BATCH/processing/log.out | sed 's/-d: //')

if [ ! -f $pth/*sequencing.csv ]; then 
  echo "ERROR: could not find lab metadata. Searched for '$pth/*sequencing.csv', but found nothing."
  exit 1
else
  labmeta=$(find $pth/*sequencing.csv)
fi

# Run .rmd script.
Rscript -e "rmarkdown::render(input='$RMD',output_file='$PWD/$BATCH/QC/$BATCH.html',knit_root_dir='$PWD',params=list(batch='$BATCH',labmeta='$labmeta'))"
