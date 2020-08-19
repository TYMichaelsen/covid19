#!/usr/bin/env bash
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-m file -s file -o dir -t int] 
-- COVID-19 pipeline for nextstrain visualization and basic genomic analysis v. $VERSION:  

Arguments:
    -h  Show this help text.
    -m  Metadata file.
    -s  Sequence file.
    -o  (Develop only) Specify output directory.
    -t  (Develop only) Number of threads.
    -f  (Develop only) Force override existing output directory. 
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hfm:s:o:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) META=$OPTARG;;
    s) SEQS=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    f) FORCE=1;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Setup directories.
THISDIR=$(dirname $(readlink -f $0))
DISTDIR="/srv/rbd/covid19"
DEPEND_DIR="${THISDIR}/dependencies"
REF="${DEPEND_DIR}/ref"


### TESTING ################################################
#THISDIR=/srv/rbd/covid19/git/covid19/workflow
#DISTDIR="/srv/rbd/tym/test-nextstrain"
#OUTDIR=testing
#DEPEND_DIR="${THISDIR}/dependencies"
#REF="${DEPEND_DIR}/ref"

#THREADS=100

#META=/srv/rbd/covid19/genomes/2020-08-17-11-59_export/metadata.tsv
#SEQS=/srv/rbd/covid19/genomes/2020-08-17-11-59_export/sequences.fasta

############################################################

# Check missing arguments
MISSING="is missing but required. Exiting."
# if [ -z ${META+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
# if [ -z ${SEQS+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then OUTDIR=${NEXTSTRAINOUT}/$(date +%Y-%m-%d)_nextstrain; fi;
# if [ -z ${OUTDIR+x} ]; then OUTDIR=$PWD/$(date +%Y-%m-%d)_nextstrain; fi;
if [ -z ${THREADS+x} ]; then THREADS=64; fi;

### Code.----------------------------------------------------------------------

# Source utility functions
source ${THISDIR}/utils.sh

# setup output folders.
#if [ -d $OUTDIR -a  x$FORCE == x  ]; then
#    echo "$OUTDIR already exists!"
#    echo "Please choose a different output directory or use -f to force delete $OUTDIR"
#    exit 1
#fi

#if [ -d $OUTDIR -a x${FORCE} != x  ]; then
#    echo "Deleting $OUTDIR ..."
#    rm -rf $OUTDIR
#fi

mkdir -p $OUTDIR
OUTDIR=$(readlink -f $OUTDIR) # Convert to absolute path
mkdir -p $OUTDIR/data

if [ -z ${META+x} ]; then
    META=$(findTheLatest "${DISTDIR}/genomes/*export")/metadata.tsv
    echo "WARNING: -m not provided, will use the latest one for nextstrain:"
    echo $META
fi
if [ -z ${SEQS+x} ]; then
    SEQS=$(findTheLatest "${DISTDIR}/genomes/*export")/sequences.fasta
    echo "WARNING: -s not provided, will use the latest sequences in :"
    echo $SEQS
fi


if [ ! -f $OUTDIR/data/metadata_nextstrain.tsv ]; then

    ###############################################################################
    # Merge the metadata with GISAID metadata.-------------------------------------
    ###############################################################################

    GISAID_META=$(findTheLatest "${DISTDIR}/global_data/*tsv")

    Rscript ${THISDIR}/merge_clean_metadata.R -l $META -g $GISAID_META -o $OUTDIR
fi
if [ ! -f $OUTDIR/data/masked.fasta ]; then

    ###############################################################################
    # Merge sequences.-------------------------------------------------------------
    ###############################################################################

    GISAID_FASTA=$(findTheLatest "${DISTDIR}/global_data/*fasta")


    # Merge GISAID and SSI sequences.
    cat $SEQS $GISAID_FASTA |
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - > $OUTDIR/masked.fasta # make one-line fasta.

    # List the sequences to include.
    awk 'NR > 1 {print $1}' $OUTDIR/metadata_nextstrain.tsv > $OUTDIR/include.txt

    # Subset sequences.
    cat $OUTDIR/masked.fasta | 
        awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
        awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' - | # tabularize.
        awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if (">"$1 in seqs) {print seqs[">"$1]} else {print $0" had no matching sequence, excluding." > "/dev/stderr"}}' - $OUTDIR/include.txt | # Subset to ones in metadata.
        tr "\t" "\n" > tmp && mv tmp $OUTDIR/masked.fasta # de-tabularize.

    # Move stuff to /data.
    mv $OUTDIR/masked.fasta $OUTDIR/data/masked.fasta
    mv $OUTDIR/metadata_full.tsv $OUTDIR/data/metadata_full.tsv
    mv $OUTDIR/metadata_nextstrain.tsv $OUTDIR/data/metadata_nextstrain.tsv
    mv $OUTDIR/include.txt $OUTDIR/data/include.txt

fi

###############################################################################
# Run nextstrain 
###############################################################################

NCOV_ROOT="/srv/rbd/bin/ncov.1308"
AUGUR_ENV="/srv/rbd/bin/conda/envs/augur"

# cd $OUTDIR
# rsync -avzp --exclude .git --exclude benmarks --exclude .snakemake --exclude .github  $NCOV_ROOT/ ./

cd $NCOV_ROOT
source activate  $AUGUR_ENV
snakemake --cores $THREADS --profile my_profiles/denmark --config metadata=$OUTDIR/data/metadata_nextstrain.tsv sequences=$OUTDIR/data/masked.fasta 

# Move results and auspice directories to $OUTDIR when snakemake finishes successfully
if [ $? -eq 0 ]; then
    cp -r results auspice $OUTDIR
    $NCOV_ROOT/scripts/assign_clades.py --nthreads $THREADS  \
                                       --alignment $OUTDIR/results/masked.fasta \
                                       --clades $OUTDIR/results/DenmarkOnly/temp_subclades.tsv \
                                       --chunk-size  $THREADS \
                                       --output $OUTDIR/results/global_clades_assignment.tsv
fi
