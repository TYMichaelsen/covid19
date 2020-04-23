#!/bin/bash

VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-m file -s file -o dir -t string] 
-- COVID-19 analysis pipeline v. $VERSION:

Arguments:
    -h  Show this help text.
    -m  File (.tsv) with metadata.
    -s  Sequences corresponding to metadata.
    -o  Output directory.
    -t  Number of threads.

Output:
    To come.

Note: 
This is beta software!
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hm:s:o:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) META=$OPTARG;;
    s) SEQS=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${META+x} ]; then echo "WARNING: -m not provided, parts of the pipeline will not be performed."; fi;
if [ -z ${SEQS+x} ]; then echo "-s $MISSING"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then OUTDIR=$PWD; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------

# setup output folders.
if [ "$OUTDIR" != "$PWD" ]; then rm -rf $OUTDIR; mkdir $OUTDIR; fi

# Capture stdout to log file.
#exec 3>&1 4>&2
#trap 'exec 2>&4 1>&3' 0 1 2 3
#exec 1>$OUTDIR/log.out 2>&1

# Settings
REF=/srv/rbd/covid19/data/reference
ROOT_SEQS=/srv/rbd/covid19/nextstrain-analysis/GISAID-data/root.fasta
ROOT_META=/srv/rbd/covid19/nextstrain-analysis/GISAID-data/root.tsv

# Add root sequences (See folder GISAID-data for details).
cat $SEQS $ROOT_SEQS > $OUTDIR/raw.fasta

### Alignment ###
# Do alignment in chunks.
rm -rf $OUTDIR/split_fasta; mkdir $OUTDIR/split_fasta
rm -rf $OUTDIR/split_align; mkdir $OUTDIR/split_align

split -d -l 40 --additional-suffix=.fasta $OUTDIR/raw.fasta $OUTDIR/split_fasta/

parallel -j$THREADS \
'
augur align \
--sequences {2}/split_fasta/{1} \
--reference-sequence {3} \
--output {2}/split_align/{1} \
--nthreads 1 \
--remove-reference \
--fill-gaps &> {2}/log.out

' ::: `ls $OUTDIR/split_fasta` ::: $OUTDIR ::: $REF/MN908947.3.gb

cat $OUTDIR/split_align/*.fasta > $OUTDIR/aligned.fasta
cat $OUTDIR/split_align/*.log > $OUTDIR/aligned.log

rm -r $OUTDIR/split_fasta
rm -r $OUTDIR/split_align

### Mask bases ###
mask_sites="18529 29849 29851 29853"

python3 /srv/rbd/covid19/nextstrain-analysis/GISAID-data/ncov/scripts/mask-alignment.py \
    --alignment $OUTDIR/aligned.fasta \
    --mask-from-beginning 130 \
    --mask-from-end 50 \
    --mask-sites $mask_sites \
    --output $OUTDIR/masked.fasta

### Haplotyping. ### 
#Rscipt haplotyping.R $OUTDIR/masked.fasta
    
### build tree. ###
augur tree \
    --alignment $OUTDIR/masked.fasta \
    --output $OUTDIR/tree_raw.nwk \
    --nthreads $THREADS

if [ -z ${META+x} ] || [ ! -f $META ]; then
  echo "No metadata provided or not founds, stopping after alignment and initial tree."
  exit 1
fi

if [ ! -f $META ]; then
  echo "ERROR: Metadata not found, exiting."; exit 1
else
  cat $META $ROOT_META > $OUTDIR/metadata.tsv
fi

### build time-adjusted tree (THE TIME KILLER!).
augur refine \
--tree $OUTDIR/tree_raw.nwk \
--alignment $OUTDIR/masked.fasta \
--metadata $OUTDIR/metadata.tsv \
--output-tree $OUTDIR/tree.nwk \
--output-node-data $OUTDIR/branch_lengths.json \
--timetree \
--root Wuhan-Hu-1/2019 Wuhan/WH01/2019 \
--clock-rate 0.0008 \
--clock-std-dev 0.0004 \
--coalescent skyline \
--date-inference marginal \
--divergence-unit mutations \
--date-confidence \
--no-covariance \
--clock-filter-iqd 4

### ancestral tree.
augur ancestral \
  --tree $OUTDIR/tree.nwk \
  --alignment $OUTDIR/masked.fasta \
  --output-node-data $OUTDIR/nt_muts.json \
  --output-sequences $OUTDIR/ancestral.fasta \
  --inference joint \
  --infer-ambiguous

### Translate NT ot AA.
augur translate \
  --tree $OUTDIR/tree.nwk \
  --ancestral-sequences $OUTDIR/nt_muts.json \
  --reference-sequence $REF/MN908947.3.gb \
  --output-node-data $OUTDIR/aa_muts.json
            
### traits.
#augur traits \
#  --tree results/tree.nwk \
#  --metadata $META \
#  --weights ncov/config/weight.tsv \
#  --output results/traits.json \
#  --columns country_exposure \
#  --confidence \
#  --sampling-bias-correction 2.5 
            
### add clades.
augur clades \
  --tree $OUTDIR/tree.nwk \
  --mutations $OUTDIR/nt_muts.json $OUTDIR/aa_muts.json \
  --clades GISAID-data/ncov/config/clades.tsv \
  --output-node-data $OUTDIR/clades.json
  
### construct colouring.
python3 GISAID-data/ncov/scripts/assign-colors.py \
  --ordering GISAID-data/ncov/config/ordering.tsv \
  --color-schemes GISAID-data/ncov/config/color_schemes.tsv \
  --output $OUTDIR/colors.tsv
  
### Construct output for auspice.
mkdir -p $OUTDIR/auspice

# Cat DK lat long with auspice.
cat GISAID-data/ncov/config/lat_longs.tsv /srv/rbd/covid19/data/metadata_SSI/latlong_DK_nextstrain.tsv > $OUTDIR/latlongs.tsv

augur export v2 \
  --tree $OUTDIR/tree.nwk \
  --metadata $OUTDIR/metadata.tsv \
  --node-data $OUTDIR/branch_lengths.json $OUTDIR/nt_muts.json $OUTDIR/aa_muts.json $OUTDIR/clades.json \
  --auspice-config GISAID-data/ncov/config/auspice_config.json \
  --colors $OUTDIR/colors.tsv \
  --color-by-metadata region country division location \
  --lat-longs $OUTDIR/latlongs.tsv \
  --output $OUTDIR/auspice/ncov_custom.json
  

#--auspice-config {input.auspice_config} \
#--description {input.description} \