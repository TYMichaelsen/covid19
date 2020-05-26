#!/bin/bash
# By Thomas Y. Michaelsen
# Heavily inspired by the ncov19 workflow provided by the nextstrain team.

set -e
# set -x

VERSION=0.2.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-m file -s file -o dir -t string] [-f (force delete)]
-- COVID-19 analysis pipeline v. $VERSION:

Arguments:
    -h  Show this help text.
    -m  File (.tsv) with metadata formatted specific for nextstrain.
    -s  Sequences corresponding to metadata.
    -o  Output directory.
    -t  Number of threads.
    -f  Force remove existing data (Be careful!)
    -c  Clades table. Will calculate if missing

Output:
    To come.

Note: 
This is beta software!
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hfm:s:o:t:c:' OPTION; do
    case $OPTION in
        h) echo "$USAGE"; exit 1;;
        m) META=$OPTARG;;
        s) SEQS=$OPTARG;;
        o) OUTDIR=$OPTARG;;
        t) THREADS=$OPTARG;;
        c) CLADES=$OPTARG;;
        f) FORCE=1;;
        :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
        \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
    esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${META+x} ]; then echo "WARNING: -m not provided, parts of the pipeline will not be performed."; fi;
if [ -z ${SEQS+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then OUTDIR=$PWD/$(date +%Y-%m-%d)_nextstrain; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------

# setup output folders.
if [ -d $OUTDIR -a  x$FORCE == x  ]; then
    echo "$OUTDIR already exists!"
    echo "Please choose a different output directory or use -f to force delete $OUTDIR"
    exit 1
fi

if [ -d $OUTDIR -a x${FORCE} != x ]; then
    echo "Deleting $OUTDIR ..."
    rm -rf $OUTDIR
fi

if [ ! -d "$OUTDIR" ]; then  mkdir $OUTDIR; fi

## Capture stdout to log file.
# exec 3>&1 4>&2
# trap 'exec 2>&4 1>&3' 0 1 2 3
# exec 1>$OUTDIR/log.out 2>&1

# Settings
# Distribution directory, the root of all input and final output
THISDIR=$(dirname $(readlink -f $0))
DEFDIR="${THISDIR}/dependencies"
DISTDIR="/srv/rbd/covid19"

REF="${DEFDIR}"
METADIR="${DISTDIR}/metadata"
NCOVDIR="${DEFDIR}/nextstrain" # Keep this for flexibility
ROOTSEQSDIR="${DEFDIR}/nextstrain"
ROOT_SEQS=${ROOTSEQSDIR}/root.fasta
ROOT_META=${ROOTSEQSDIR}/root.tsv

# Add root sequences (See folder GISAID-data for details).
cat $SEQS $ROOT_SEQS > $OUTDIR/raw.fasta

# Append \t to root meta to match metadata.
rootcols=$(awk 'NR == 1 {print NF}' $META)
cat $META <(awk -F'\t' -v Nroot=$rootcols 'BEGIN {OFS="\t"} { for(i=NF+1; i<=Nroot; i++) $i=""; print}' $ROOT_META) > $OUTDIR/metadata.tsv
#cat $META $ROOT_META > $OUTDIR/metadata.tsv

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

if [ -z ${META+x} ] || [ ! -f $META ]; then
    echo "No metadata provided or not founds, stopping after alignment and initial tree."
    exit 1
fi

if [ ! -f $META ]; then
    echo "ERROR: Metadata not found, exiting."; exit 1
else
    cat $META $ROOT_META > $OUTDIR/metadata.tsv
    # cat $META > $OUTDIR/metadata.tsv
    # Replace fasta header (library_id) withs strain names (ssi_id)
    # awk -F'\t' '
    # (FNR==NR){
    #     lib2id[$3]=$8; next}
    # { if ($0 ~/^>/) { hdr=$0; sub(">","", hdr);
    #         if (length(lib2id[hdr]) >0)
    #         {$0=">"lib2id[hdr]; }
    #  }
    #  print $0} ' \
        # ${METADIR}/2020-04-28-19-57_metadata.tsv $OUTDIR/aligned.fasta > ${OUTDIR}/aligned.fixheader.fasta
fi

### Mask bases ###
mask_sites="18529 29849 29851 29853"

python3 ${DEFDIR}/mask-alignment.py \
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

### Make clades.tsv from data of Lineage repo 
if [ ! -f "${CLADES}" ]; then
    bash -c "source activate pangolin; python3 ${DEFDIR}/make_clade_tsv.py --clade-out=${OUTDIR}/clades.tsv"
else
    cp ${CLADES} ${OUTDIR}/clades.tsv
fi

### add clades.
augur clades \
      --tree $OUTDIR/tree.nwk \
      --mutations $OUTDIR/nt_muts.json $OUTDIR/aa_muts.json \
      --clades ${OUTDIR}/clades.tsv \
      --output-node-data $OUTDIR/clades.json

### construct colouring.
python3 ${NCOVDIR}/assign-colors.py \
        --ordering ${NCOVDIR}/ordering.tsv \
        --color-schemes ${NCOVDIR}/color_schemes.tsv \
        --output $OUTDIR/colors.tsv \
        --metadata ${OUTDIR}/metadata.tsv

### Construct frequency tables.
augur frequencies \
      --method kde \
      --metadata $OUTDIR/metadata.tsv \
      --tree $OUTDIR/tree.nwk \
      --alignments $OUTDIR/masked.fasta \
      --output $OUTDIR/tip-frequencies.json

### Running pangolin
# Pangolin itself has a conda environment.
# And we are running inside augur conda one.
# Therefore it is neccessary to run pangolin
# this way (wrap inside a `bash -c` session)

# Pangolin creates quite large temporary files
### Doesn't work at the moment.
# See this issue: https://github.com/hCoV-2019/pangolin/issues/61
bash -c "
source activate pangolin
pangolin $OUTDIR/masked.fasta -t $THREADS \
    --tempdir pangtmp \
    --outdir $OUTDIR \
    --include-putative
    rm -rf ./pangtmp
"

# Add pangolin lineages and lineage version to metadata.
join -1 1 -2 1 -t $'\t' <(tail -n +2 $OUTDIR/metadata.tsv | sort) <(awk -F',' '{print $1"\t"$2"\t"$5}' $OUTDIR/lineage_report.csv | tail -n +2 | sort) > $OUTDIR/metadata_w_linage.tsv

cat <(awk 'NR == 1 {print $0"\tlineage\tlineage_version"}' $OUTDIR/metadata.tsv) $OUTDIR/metadata_w_linage.tsv > tmp && mv tmp $OUTDIR/metadata_w_linage.tsv

### Export data for auspice
# Construct output for auspice.
mkdir -p $OUTDIR/auspice

# Cat DK lat long with auspice.
cat ${METADIR}/latlong_nextstrain.tsv ${NCOVDIR}/lat_longs.tsv > $OUTDIR/latlongs.tsv

# Get all metadata columns except strain, virus and date to display in auspice.
cols=$(awk -F'\t' 'NR == 1 {$1=$2=$3=""; print $0}' $OUTDIR/metadata_w_linage.tsv)

augur export v2 \
      --tree $OUTDIR/tree.nwk \
      --metadata $OUTDIR/metadata_w_linage.tsv \
      --node-data $OUTDIR/branch_lengths.json $OUTDIR/nt_muts.json $OUTDIR/aa_muts.json $OUTDIR/clades.json \
      --auspice-config ${NCOVDIR}/auspice_config_AAU.json \
      --color-by-metadata $cols \
      --lat-longs $OUTDIR/latlongs.tsv \
      --output $OUTDIR/auspice/ncov_custom.json
#--description {input.description} \

cp $OUTDIR/tip-frequencies.json $OUTDIR/auspice/ncov_custom_tip-frequencies.json

