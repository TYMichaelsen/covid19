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
    -p  (Develop only) Pull ncov from github
    -k  (Develop only) Additional snakemake arguments given in quotes

"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hfpm:s:o:t:k:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) META=$OPTARG;;
    s) SEQS=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    f) FORCE=1;;
    p) PULL_GITHUB=1;;
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
# if [ -z ${META+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
# if [ -z ${SEQS+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then OUTDIR=${NEXTSTRAINOUT}/$(date +%Y-%m-%d)_nextstrain; fi;
# if [ -z ${OUTDIR+x} ]; then OUTDIR=$PWD/$(date +%Y-%m-%d)_nextstrain; fi;
if [ -z ${THREADS+x} ]; then THREADS=64; fi;

### Code.----------------------------------------------------------------------
# Source utility functions
source ${THISDIR}/utils.sh

IMGDIR="/srv/rbd/thecontainer"
SINGIMG=$(findTheLatest "${IMGDIR}/*sif")
# SINGIMG=singularity/covid19_latest.sif
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
ARGSTR="--cores $THREADS --profile my_profiles/denmark --config metadata=$OUTDIR/data/metadata_nextstrain.tsv sequences=$OUTDIR/data/masked.fasta ${SNAKE_ADD}"

# Run nextstrain 
###############################################################################

if [ -d $OUTDIR -a x${PULL_GITHUB} != x  ]; then
    cd $OUTDIR
    wget https://github.com/biocyberman/ncov/archive/aau.zip
    unzip aau.zip && rm aau.zip
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

    # Move results and auspice directories to $OUTDIR when snakemake finishes successfully
    if [ $? -eq 0 ]; then
        cp -r results auspice $OUTDIR
        $NCOV_ROOT/scripts/assign_clades.py --nthreads $THREADS  \
                                            python $NCOV_ROOT/scripts/assign_clades.py --nthreads $THREADS  \
                                            --alignment $OUTDIR/results/masked.fasta \
                                            --clades $OUTDIR/results/DenmarkOnly/temp_subclades.tsv \
                                            --chunk-size  $THREADS \
                                            --output $OUTDIR/results/global_clades_assignment.tsv
    fi
else
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

# Move results and auspice directories to $OUTDIR when snakemake finishes successfully
if [ -f auspice/ncov_DenmarkGlobal.json  ]; then
    echo Copying output to $OUTDIR ...
    cp -r results auspice $OUTDIR
    echo Running assign_clades ...
    python $NCOV_ROOT/scripts/assign_clades.py --nthreads $THREADS  \
                                       --alignment $OUTDIR/results/masked.fasta \
                                       --clades $OUTDIR/results/DenmarkOnly/temp_subclades.tsv \
                                       --chunk-size  $THREADS \
                                       --output $OUTDIR/results/global_clades_assignment.tsv
fi
HEREDOC

fi
