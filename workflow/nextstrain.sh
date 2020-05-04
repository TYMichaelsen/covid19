#!/bin/bash
# By Thomas Y. Michaelsen
# Heavily inspired by the ncov19 workflow provided by the nextstrain team.

set -ex

VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-m file -s file -o dir -t string] 
-- COVID-19 analysis pipeline v. $VERSION:

Arguments:
    -h  Show this help text.
    -m  File (.tsv) with metadata formatted specific for nextstrain.
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
if [ -z ${SEQS+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
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
# Distribution directory, the root of all input and final output
DISTDIR="/srv/rbd/covid19/current"

REF="${DISTDIR}/auxdata/reference/"
METADIR="${DISTDIR}/metadata/metadata_SSI"
NCOVDIR="${DISTDIR}/auxdata/ncov"
ROOTSEQSDIR="${DISTDIR}/auxdata/root_seqs"
ROOT_SEQS=${ROOTSEQSDIR}/root.fasta
ROOT_META=${ROOTSEQSDIR}/root.tsv

# Add root sequences (See folder GISAID-data for details).
cat $SEQS $ROOT_SEQS > $OUTDIR/raw.fasta
cat $META $ROOT_META > $OUTDIR/metadata.tsv

### Alignment ###
# Do alignment in chunks.
rm -rf $OUTDIR/split_fasta; mkdir $OUTDIR/split_fasta
rm -rf $OUTDIR/split_align; mkdir $OUTDIR/split_align

split -d -l 40 --additional-suffix=.fasta $OUTDIR/raw.fasta $OUTDIR/split_fasta/
FASTAFILES=`ls $OUTDIR/split_fasta` 

parallel -j$THREADS \
'
augur align \
--sequences {2}/split_fasta/{1} \
--reference-sequence {3} \
--output {2}/split_align/{1}.aligned \
--nthreads 1 \
--remove-reference \
--fill-gaps &> {2}/log.out

' ::: $FASTAFILES ::: $OUTDIR ::: $REF/MN908947.3.gb

cat $OUTDIR/split_align/*.aligned > $OUTDIR/aligned.fasta
cat $OUTDIR/split_align/*.log > $OUTDIR/aligned.log

rm -r $OUTDIR/split_fasta
rm -r $OUTDIR/split_align

### Mask bases ###
mask_sites="18529 29849 29851 29853"

python3 ${NCOVDIR}/scripts/mask-alignment.py \
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
augur traits \
 --tree ${OUTDIR}/tree.nwk \
 --metadata ${OUTDIR}/metadata.tsv \
 --output ${OUTDIR}/traits.json \
 --columns country_exposure \
 --confidence \
 --sampling-bias-correction 2.5 
            
### add clades.
augur clades \
  --tree $OUTDIR/tree.nwk \
  --mutations $OUTDIR/nt_muts.json $OUTDIR/aa_muts.json \
  --clades ${ROOTSEQSDIR}/pangolin_clades.tsv \
  --output-node-data $OUTDIR/clades.json
  
### construct colouring.
python3 ${NCOVDIR}/scripts/assign-colors.py \
  --ordering ${NCOVDIR}/config/ordering.tsv \
  --color-schemes ${NCOVDIR}/config/color_schemes.tsv \
  --output $OUTDIR/colors.tsv

### Construct frequency tables.
augur frequencies \
--method kde \
--metadata $OUTDIR/metadata.tsv \
--tree $OUTDIR/tree.nwk \
--alignments $OUTDIR/masked.fasta \
--output $OUTDIR/tip-frequencies.json

### Construct output for auspice.
mkdir -p $OUTDIR/auspice

# Cat DK lat long with auspice.
cat ${NCOVDIR}/config/lat_longs.tsv ${METADIR}/latlong_DK_nextstrain.tsv > $OUTDIR/latlongs.tsv

augur export v2 \
  --tree $OUTDIR/tree.nwk \
  --metadata $OUTDIR/metadata.tsv \
  --node-data $OUTDIR/branch_lengths.json $OUTDIR/nt_muts.json $OUTDIR/aa_muts.json ${OUTDIR}/traits.json $OUTDIR/clades.json $OUTDIR/tip-frequencies.json \
  --auspice-config ${NCOVDIR}/config/auspice_config.json \
  --colors $OUTDIR/colors.tsv \
  --color-by-metadata region country division location \
  --lat-longs $OUTDIR/latlongs.tsv \
  --panels tree map entropy frequencies \
  --output $OUTDIR/auspice/ncov_custom.json
  

#--auspice-config {input.auspice_config} \
#--description {input.description} \
