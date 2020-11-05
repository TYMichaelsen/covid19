#!/usr/bin/env bash
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-m file -s file -o dir -t int]
-- Exstract SNPS directly from raw input fasta file v. $VERSION:  

Arguments:
    -h  Show this help text.
    -s  Sequence file.
    -b  Choose custom build among 'light' (default), 'full' or 'global'. Multiple builds also possible: light,full (no space around ',')
    -i  (Develop only) Specify a specific singularity image to use.
    -o  (Develop only) Specify output directory.
    -k  (Develop only) Additional snakemake arguments given in quotes

"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hfpm:b:s:o:t:k:i:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    s) SEQS=$OPTARG;;
    m) META=$OPTARG;;
    b) CUSTOM_BUILD=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    i) SINGIMG=$OPTARG;;
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

### Code.----------------------------------------------------------------------
# Source utility functions
# source ${THISDIR}/utils.sh
source ${DISTDIR}/git/covid19/workflow/utils.sh

IMGDIR="/srv/rbd/thecontainer"
if [ -z "$SINGIMG" ]; then
    SINGIMG=$(findTheLatest "${IMGDIR}/*sif")
fi
echo Using Singularity image: $SINGIMG

if [ -z ${SEQS+x} ]; then
    SEQS=$(findTheLatest "${DISTDIR}/genomes/*export")/sequences.fasta
    echo "WARNING: -s not provided, will use the latest sequences in :"
    echo $SEQS
fi

if [ -z ${META+x} ]; then
    META=$(findTheLatest "${DISTDIR}/metadata/*nextstrain.tsv")
    echo "WARNING: -m not provided, will use the latest one for nextstrain:"
    echo $META
fi

if [ -z ${OUTDIR+x} ]; then 
    OUTDIR=${NEXTSTRAINOUT}/$(basename -s _export $(findTheLatest "${DISTDIR}/genomes/*export"))_nextstrain/mutations
fi

if [ -z ${CUSTOM_BUILD+x} ]; then 
    CUSTOM_BUILD="light"
fi

if [ -n "$CUSTOM_BUILD" ]; then
    SNAKE_ADD="custom_build=$CUSTOM_BUILD $SNAKE_ADD"
fi



GENOMEDIR=$(dirname $SEQS)


mkdir -p $OUTDIR
OUTDIR=$(readlink -f $OUTDIR) # Convert to absolute path
mkdir -p $OUTDIR/data
RUNLOGFILE=${OUTDIR}/'mutations_run.log'

GISAID_META=$(findTheLatest "${DISTDIR}/global_data/*tsv")
GISAID_FASTA=$(findTheLatest "${DISTDIR}/global_data/*fasta")
NCOV_ROOT="/opt/nextstrain/ncov-aau"
BPROFILE="my_profiles/denmark"

gisaid_date=$(basename $GISAID_META|sed -e 's/metadata_//;s/.tsv//')
data_date=$(basename -s "_export" $(dirname $SEQS))


# Make sure when re-running on the same $OUTDIR use the same input data
SNAKE_OUT=${OUTDIR}
if [ -s "${RUNLOGFILE}" ]; then
    data_date=$(grep snakemake ${RUNLOGFILE} |tail -n 1|sed 's/ /\n/g'| grep 'data_date=' | awk -F '=' '{print $2}'|tr -d \'\")
    gisaid_date=$(grep snakemake ${RUNLOGFILE} |tail -n 1|sed 's/ /\n/g'| grep 'gisaid_date=' | awk -F '=' '{print $2}'|tr -d \'\")

    META=$(grep snakemake ${RUNLOGFILE} |tail -n 1|sed 's/ /\n/g'| grep 'denmark_meta=' | awk -F '=' '{print $2}')
    SEQS=$(grep snakemake ${RUNLOGFILE} |tail -n 1|sed 's/ /\n/g'| grep 'denmark_fasta=' | awk -F '=' '{print $2}')
    GISAID_META=$(grep snakemake ${RUNLOGFILE} |tail -n 1|sed 's/ /\n/g'| grep 'gisaid_meta=' | awk -F '=' '{print $2}')
    GISAID_FASTA=$(grep snakemake ${RUNLOGFILE} |tail -n 1|sed 's/ /\n/g'| grep 'gisaid_fasta=' | awk -F '=' '{print $2}')
    echo Existing run output detected at $OUTDIR.
    echo Now rerunning with data_date: $data_date, gisaid_date: $gisaid_date ...
fi

THREADS=48
SMKFILE="${BPROFILE}/extract_SNPs.smk"
ARGSTR="--snakefile=${SMKFILE} --cores $THREADS --profile $BPROFILE "
ARGSTR="$ARGSTR --config gisaid_date='\"$gisaid_date\"' data_date='\"$data_date\"' "
ARGSTR="$ARGSTR denmark_meta=${META} denmark_fasta=${SEQS} gisaid_meta=${GISAID_META} gisaid_fasta=${GISAID_FASTA}"
ARGSTR="$ARGSTR outdir=${SNAKE_OUT}  ${SNAKE_ADD}"

# Run extract mutations
###############################################################################
TimeStamp=$(date +%y%m%d_%H%M)

echo ${TimeStamp}: ${USER}, PWD: ${PWD} >>${RUNLOGFILE}
echo ${TimeStamp}: "$0 $@" >>${RUNLOGFILE}
echo ${TimeStamp}: Singularity Image Used: $SINGIMG >> ${RUNLOGFILE}
echo ${TimeStamp}: snakemake ${ARGSTR} >> ${RUNLOGFILE}
singularity exec  -B /srv/rbd:/srv/rbd \
            -B $HOME:$HOME \
            $SINGIMG bash <<HEREDOC

source activate nextstrain

###############################################################################
# Run extract mutations 
###############################################################################
if [ ! -d $OUTDIR/ncov-aau ]; then
   cp -r $NCOV_ROOT $OUTDIR
fi
cd $OUTDIR/ncov-aau
snakemake ${ARGSTR}
HEREDOC
