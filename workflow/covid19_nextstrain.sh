#!/usr/bin/env bash
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-m file -s file -o dir -t int]
-- COVID-19 pipeline for nextstrain visualization and basic genomic analysis v. $VERSION:  

Arguments:
    -h  Show this help text.
    -m  Metadata file.
    -s  Sequence file.
    -g  Flag to build global dataset. Default is False, not to build, only take all Danish data and subsample global data
    -i  (Develop only) Specify a specific singularity image to use.
    -o  (Develop only) Specify output directory.
    -t  (Develop only) Number of threads.
    -f  (Develop only) Force override existing output directory. 
    -p  (Develop only) Pull ncov from github
    -k  (Develop only) Additional snakemake arguments given in quotes

"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hfpgm:s:o:t:k:i:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) META=$OPTARG;;
    s) SEQS=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    i) SINGIMG=$OPTARG;;
    f) FORCE=1;;
    p) PULL_GITHUB=1;;
    g) BUILD_GRLOBAL=1;;
    k) SNAKE_ADD=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Setup directories.
THISDIR=$(dirname $(readlink -f $0))
DISTDIR="/srv/rbd/covid19"
NEXTSTRAINOUT="${DISTDIR}/nextstrain"

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${THREADS+x} ]; then THREADS=64; fi;

### Code.----------------------------------------------------------------------
# Source utility functions
source ${THISDIR}/utils.sh

IMGDIR="/srv/rbd/thecontainer"
if [ -z "$SINGIMG" ]; then
    SINGIMG=$(findTheLatest "${IMGDIR}/*sif")
fi
echo Using Singularity image: $SINGIMG
# DISTDIR="/srv/rbd/covid19/current"

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

if [ -z ${OUTDIR+x} ]; then 
  OUTDIR=$(basename -s _export $(findTheLatest "${DISTDIR}/genomes/*export"))_nextstrain
fi

GENOMEDIR=$(dirname $SEQS)
METADIR=$(dirname $META)

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

GISAID_META=$(findTheLatest "${DISTDIR}/global_data/*tsv")
GISAID_FASTA=$(findTheLatest "${DISTDIR}/global_data/*fasta")
NCOV_ROOT="/opt/nextstrain/ncov-aau"
if [ -n "$BUILD_GLOBAL" ]; then
    BPROFILE="my_profiles/denmark"
else
    BPROFILE="my_profiles/denmarkonly"
fi


gisaid_date=$(basename $GISAID_META|sed -e 's/metadata_//;s/.tsv//')
data_date=$(basename -s "_export" $(dirname $SEQS))

ARGSTR="--cores $THREADS --profile $BPROFILE --config gisaid_date=$gisaid_date data_date=$data_date "
ARGSTR="$ARGSTR metadata=$OUTDIR/data/metadata_nextstrain.tsv sequences=$OUTDIR/data/masked.fasta "
ARGSTR="$ARGSTR outdir=${OUTDIR}/results out_auspice=${OUTDIR}/aupsice ${SNAKE_ADD}"

# Run nextstrain 
###############################################################################

if [ -d $OUTDIR -a x${PULL_GITHUB} != x ]; then
    cd $OUTDIR
    wget https://github.com/biocyberman/ncov/archive/aau.zip
    unzip -o aau.zip && rm aau.zip
    cd -
fi

CONDA_RUN=0 # preserve conda-based setup for convenience. 

if [ $CONDA_RUN -eq 1 ]; then

    NCOV_ROOT="/srv/rbd/bin/ncov.1308"
    AUGUR_ENV="/srv/rbd/bin/conda/envs/augur"

    # cd $OUTDIR
    # rsync -avzp --exclude .git --exclude benmarks --exclude .snakemake --exclude .github  $NCOV_ROOT/ ./
    if [ ! -f $OUTDIR/data/metadata_nextstrain.tsv ]; then
        ###############################################################################
        # Merge the metadata with GISAID metadata.-------------------------------------
        ###############################################################################

        Rscript --vanilla --no-environ ${THISDIR}/merge_clean_metadata.R -l $META -g $GISAID_META -o $OUTDIR
        mv $OUTDIR/metadata_nextstrain.tsv $OUTDIR/data/metadata_nextstrain.tsv
        mv $OUTDIR/metadata_full.tsv $OUTDIR/data/metadata_full.tsv
    fi
    if [ ! -f $OUTDIR/data/masked.fasta ]; then
        ###############################################################################
        # Merge sequences.-------------------------------------------------------------
        ###############################################################################

        # Merge GISAID and SSI sequences.
        # List the sequences to include.
        awk 'NR > 1 {print \$1}' $OUTDIR/data/metadata_nextstrain.tsv > $OUTDIR/include.txt
        cat $SEQS $GISAID_FASTA | seqtk subseq  - $OUTDIR/include.txt > $OUTDIR/masked.fasta

        # Move stuff to /data.
        mv $OUTDIR/masked.fasta $OUTDIR/data/masked.fasta
        mv $OUTDIR/include.txt $OUTDIR/data/include.txt
    fi

    cd $NCOV_ROOT
    source activate  $AUGUR_ENV
    snakemake ${ARGSTR}

else
    echo Singularity Image: $SINGIMG > ${OUTDIR}/run.log
    echo snakemake ${ARGSTR} >> $OUTDIR/run.log
    singularity exec  -B /srv/rbd:/srv/rbd \
                -B $HOME:$HOME \
                $SINGIMG bash <<HEREDOC
source activate nextstrain

# Source utility functions
source ${THISDIR}/utils.sh

if [ ! -f $OUTDIR/data/metadata_nextstrain.tsv ]; then
    ###############################################################################
    # Merge the metadata with GISAID metadata.-------------------------------------
    ###############################################################################
    Rscript --vanilla --no-environ ${THISDIR}/merge_clean_metadata.R -l $META -g $GISAID_META -o $OUTDIR
    mv $OUTDIR/metadata_nextstrain.tsv $OUTDIR/data/metadata_nextstrain.tsv
    mv $OUTDIR/metadata_full.tsv $OUTDIR/data/metadata_full.tsv
fi
if [ ! -f $OUTDIR/data/masked.fasta ]; then
    ###############################################################################
    # Merge sequences.-------------------------------------------------------------
    ###############################################################################

    # Merge GISAID and SSI sequences.
    # List the sequences to include.
    awk 'NR > 1 {print \$1}' $OUTDIR/data/metadata_nextstrain.tsv > $OUTDIR/include.txt
    cat $SEQS $GISAID_FASTA | seqtk subseq - $OUTDIR/include.txt > $OUTDIR/masked.fasta

    # Move stuff to /data.
    mv $OUTDIR/masked.fasta $OUTDIR/data/masked.fasta
    mv $OUTDIR/include.txt $OUTDIR/data/include.txt

fi
###############################################################################
# Run nextstrain 
###############################################################################
if [ ! -d $OUTDIR/ncov-aau ]; then
   cp -r $NCOV_ROOT $OUTDIR
fi
cd $OUTDIR/ncov-aau
snakemake ${ARGSTR}
HEREDOC

fi
