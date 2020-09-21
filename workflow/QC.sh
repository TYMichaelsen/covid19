#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d dir -r file] 
-- COVID-19 QC script v. $VERSION:  

Arguments:
    -h  Show this help text.
    -i  Input directory.
    -b  What batch to do QC for.
    -s  What scheme you are running. See Schemes below.
    -r  R script to run to generate QC report.
    -t  Number of threads.

Schemes:
    aau_long_v3.1
    aau_short_v3
    v1
    v2
    v3
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hi:b:s:r:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    i) INPUT_DIR=$OPTARG;;
    b) BATCH=$OPTARG;;
    s) SCHEME=$OPTARG;;
    r) RMD=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${INPUT_DIR+x} ]; then echo "-i $MISSING"; exit 1; fi;
if [ -z ${BATCH+x} ]; then echo "-b $MISSING"; exit 1; fi;
if [ -z ${SCHEME+x} ]; then echo "-s $MISSING"; exit 1; fi;
if [ -z ${RMD+x} ]; then echo "-r $MISSING"; exit 1; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------

AAU_COVID19_PATH="$(dirname "$(readlink -f "$0")")"

# Setup dirs.
mkdir -p $INPUT_DIR/QC
mkdir -p $INPUT_DIR/QC/aligntree

# Logging
LOG_NAME="$INPUT_DIR/QC/QC_log_$(date +"%Y-%m-%d_%H-%M").txt"
echo "QC log" >> $LOG_NAME
#echo "AAU COVID-19 revision - $(git -C $AAU_COVID19_PATH rev-parse --short HEAD)" >> $LOG_NAME
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

REF=$AAU_COVID19_PATH/dependencies/ref/MN908947.3.gb
CLADES=$AAU_COVID19_PATH/dependencies/nextstrain/pangolin_clades.tsv

###############################################################################
# Setup data to be used in QC.
###############################################################################

# Copy over sequences.
cp $INPUT_DIR/processing/results/consensus.fasta $INPUT_DIR/QC/aligntree/sequences.fasta
#cat export/*_export/sequences.fasta > QC/aligntree/sequences.fasta

### Alignment ###
augur align \
--sequences $INPUT_DIR/QC/aligntree/sequences.fasta \
--reference-sequence $REF \
--output $INPUT_DIR/QC/aligntree/aligned.fasta \
--nthreads $THREADS &> $INPUT_DIR/QC/aligntree/log.out

### Mask bases ###
mask_sites="18529 29849 29851 29853"

python3 $AAU_COVID19_PATH/mask-alignment.py \
--alignment $INPUT_DIR/QC/aligntree/aligned.fasta \
--mask-from-beginning 130 \
--mask-from-end 50 \
--mask-sites $mask_sites \
--output $INPUT_DIR/QC/aligntree/masked.fasta
    
### Tree ###
augur tree \
--alignment $INPUT_DIR/QC/aligntree/masked.fasta \
--output $INPUT_DIR/QC/aligntree/tree_raw.nwk \
--nthreads $THREADS
  
#augur refine \
#--tree $INPUT_DIR/QC/aligntree/tree_raw.nwk \
#--output-tree $INPUT_DIR/QC/aligntree/tree.nwk

### ancestral tree.
#augur ancestral \
#  --tree $INPUT_DIR/QC/aligntree/tree.nwk \
#  --alignment $INPUT_DIR/QC/aligntree/masked.fasta \
#  --output-node-data $INPUT_DIR/QC/aligntree/nt_muts.json \
#  --inference joint \
#  --infer-ambiguous
  
### Translate NT ot AA.
#augur translate \
#  --tree $INPUT_DIR/QC/aligntree/tree.nwk \
#  --ancestral-sequences $INPUT_DIR/QC/aligntree/nt_muts.json \
#  --reference-sequence $REF \
#  --output-node-data $INPUT_DIR/QC/aligntree/aa_muts.json
                       
### add clades.
#augur clades \
#  --tree $INPUT_DIR/QC/aligntree/tree.nwk \
#  --mutations $INPUT_DIR/QC/aligntree/nt_muts.json $INPUT_DIR/QC/aligntree/aa_muts.json \
#  --clades $CLADES \
#  --output-node-data $INPUT_DIR/QC/aligntree/clades.json
  
###############################################################################
# Generate the QC report.
###############################################################################

# Fetch path lab metadata.            
pth=$(grep "\-d" $INPUT_DIR/processing/log.out | sed 's/-d: //')

if [ ! -f $pth/sample_sheet.csv ]; then 
  echo "ERROR: could not find lab metadata. Searched for '$pth/*sequencing.csv', but found nothing."
  exit 1
else
  labmeta=$(find $pth/sample_sheet.csv)
fi

# Run .rmd script.
REF_PATH=$AAU_COVID19_PATH/dependencies/ref/MN908947.3.fasta
SCHEME_PATH=$AAU_COVID19_PATH/dependencies/primer_schemes/nCoV-2019/$SCHEME

Rscript \
  -e \
  "
  rmarkdown::render(
    input='$RMD',
    output_file='$INPUT_DIR/QC/$BATCH.html',
    knit_root_dir='$INPUT_DIR',
    params=list(
      batch='$BATCH',
      labmeta='$labmeta',
      input_dir='$INPUT_DIR',
      scheme_dir='$SCHEME_PATH',
      ref='$REF_PATH'
    )
  )
  "
