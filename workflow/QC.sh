#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.2.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d dir -r file] 
-- COVID-19 QC script v. $VERSION:  

Arguments:
    -h  Show this help text.
    -i  Input directory.
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
while getopts ':hi:s:r:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    i) INPUT_DIR=$OPTARG;;
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
echo "Command: $0 $*" >> $LOG_NAME
exec 1>>$LOG_NAME
exec 2>&1

REF=$AAU_COVID19_PATH/dependencies/ref/MN908947.3.gb

###############################################################################
# Setup data to be used in QC.
###############################################################################

# Remove failed genomes from downstream tree building and clade assignment.
MAXN=3000
MINLENGTH=25000

rm -f $INPUT_DIR/QC/aligntree/failed.fasta

cat $INPUT_DIR/processing/results/consensus.fasta |
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
awk -v THR=$MAXN -v LEN=$MINLENGTH -v outdir=$INPUT_DIR/QC/aligntree '!/^>/ { next } { getline seq; seq2=seq; Nn=gsub(/N/,"",seq) }; {if (length(seq2) > LEN && Nn <= THR) { print $0 "\n" seq2 } else {print $0 "\n" seq2 >> outdir"/failed.fasta"}}' - > $INPUT_DIR/QC/aligntree/sequences.fasta # Tidy header.

if [ -s $INPUT_DIR/QC/aligntree/failed.fasta ]; then nfail=$(grep ">" -c $INPUT_DIR/QC/aligntree/failed.fasta); else nfail=0; fi
echo "$nfail genomes were not used for tree building and nextclade, see $INPUT_DIR/QC/failed.fasta" 

source activate nextstrain

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
 
source deactivate nextstrain

### Nextclade ###

source activate nextclade 

nextclade \
--input-pcr-primers ${AAU_COVID19_PATH}/dependencies/primer_schemes/nCoV-2019/${SCHEME}/custom_primer.csv \
--input-fasta $INPUT_DIR/QC/aligntree/sequences.fasta \
--output-tsv $INPUT_DIR/QC/nextclade.tsv \
--jobs 5
 
source deactivate nextclade

### Pangolin ###

source activate pangolin 

pangolin $INPUT_DIR/QC/aligntree/sequences.fasta -t 1 --outfile $INPUT_DIR/QC/pangolin.csv --tempdir $INPUT_DIR/QC

source deactivate

###############################################################################
# Generate the QC report.
###############################################################################

if [ ! -f $INPUT_DIR/sample_sheet.csv ]; then 
  echo "ERROR: could not find lab metadata. Searched for '$INPUT_DIR/sample_sheet.csv', but found nothing."
  exit 1
else
  labmeta=$INPUT_DIR/sample_sheet.csv
fi

# Run .rmd script.
REF_PATH=$AAU_COVID19_PATH/dependencies/ref/MN908947.3.fasta
SCHEME_PATH=$AAU_COVID19_PATH/dependencies/primer_schemes/nCoV-2019/$SCHEME
BATCH=$(basename $INPUT_DIR)

Rscript \
  -e \
  "
  rmarkdown::render(
    input='$RMD',
    output_file='$INPUT_DIR/QC/$BATCH.html',
    knit_root_dir='$INPUT_DIR',
    params=list(
      labmeta='$labmeta',
      input_dir='$INPUT_DIR',
      scheme_dir='$SCHEME_PATH',
      ref='$REF_PATH'
    )
  )
  " 
