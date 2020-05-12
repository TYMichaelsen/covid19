#!/usr/bin/env bash
VERSION=0.2

# Parameter handling mirror almost exact that of nextstrain.sh
USAGE="$(basename "$0") [-h] [-m file -s file -o dir -t string] [-f (force delete)]
-- Singularity wrapper for COVID-19 analysis pipeline v. $VERSION:

Arguments:
    -h  Show this help text.
    -m  File (.tsv) with metadata formatted specific for nextstrain (Default: the latest metadata for nextstrain)
    -s  Sequences corresponding to metadata.
    -o  Output directory.
    -t  Number of threads.
    -c  Clades table. Will calculate if missing
    -f  Force remove existing data (Be careful!)

Output:
    Output of nextstrain.sh 
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hfm:s:o:t:c:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) METAFILE=$OPTARG;;
    s) GENOMEFASTA=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    c) CLADES=$OPTARG;;
    f) FORCE=1 ;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

THISDIR=$(dirname $(readlink -f $0))
DISTDIR="/srv/rbd/covid19"
METADIR="${DISTDIR}/metadata"
NEXTSTRAINOUT="${DISTDIR}/nextstrain"

# Source utility functions
source ${THISDIR}/utils.sh

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${METAFILE+x} ]; then
    METAFILE=$(findTheLatest "${METADIR}/*metadata_nextstrain.tsv")
    echo "WARNING: -m not provided, will use the latest one for nextstrain:"
    echo $METAFILE
fi
if [ -z ${GENOMEFASTA+x} ]; then
    GENOMEFASTA=$(findTheLatest "${DISTDIR}/genomes/*export")/sequences.fasta
    echo "WARNING: -s not provided, will use the latest sequences in :"
    echo $GENOMEFASTA
fi
if [ -z ${OUTDIR+x} ]; then OUTDIR=${NEXTSTRAINOUT}/$(date +%Y-%m-%d)_nextstrain; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------
DISTDIR="/srv/rbd/covid19"
METADIR="${DISTDIR}/metadata"

# setup output folders.
if [ -d $OUTDIR -a  x$FORCE == x  ]; then
    echo "$OUTDIR already exists!"
    echo "Please choose a different output directory or use -f to force delete $OUTDIR"
    exit 1
fi

if [ -d $OUTDIR -a x${FORCE} != x  ]; then
    echo "Deleting $OUTDIR ..."
    rm -rf $OUTDIR
fi


IMGDIR="/srv/rbd/thecontainer"
SINGIMG=$(findTheLatest "${IMGDIR}/*sif")
echo Using Singularity image: $SINGIMG
# DISTDIR="/srv/rbd/covid19/current"
DISTDIR="/srv/rbd"

GENOMEDIR=$(dirname $GENOMEFASTA)
METADIR=$(dirname $METAFILE)
NEXTSTRAIN_SCRIPT=${THISDIR}/nextstrain.sh

# Uncomment command below to run under augur conda envrinment for testing
# $NEXTSTRAIN_SCRIPT -s $GENOMEFASTA -m $METAFILE -o $OUTDIR -t $THREADS
# exit 1

echo -e "Running comand:\n---"
echo "singularity exec  -B $DISTDIR:$DISTDIR
-B $HOME:$HOME
-B $METADIR:$METADIR
-B $GENOMEDIR:$GENOMEDIR
$SINGIMG bash -c \"source activate nextstrain; $NEXTSTRAIN_SCRIPT -s $GENOMEFASTA -m $METAFILE -o $OUTDIR -t $THREADS\""
echo "---"


singularity exec  -B $DISTDIR:$DISTDIR \
            -B $HOME:$HOME \
            -B $METADIR:$METADIR \
            -B $GENOMEDIR:$GENOMEDIR \
            $SINGIMG bash -c "source activate nextstrain; $NEXTSTRAIN_SCRIPT -s $GENOMEFASTA -m $METAFILE -o $OUTDIR -t $THREADS"
