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
THISDIR=/srv/rbd/covid19/git/covid19/workflow
DISTDIR="/srv/rbd/tym/test-nextstrain"
OUTDIR=testing

THREADS=100

META=/srv/rbd/tym/test-nextstrain/genomes/2020-07-01-18-35_export/metadata.tsv
SEQS=/srv/rbd/tym/test-nextstrain/genomes/2020-07-01-18-35_export/sequences.fasta

############################################################

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${META+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${SEQS+x} ]; then echo "-s $MISSING"; echo "$USAGE"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then OUTDIR=$PWD/$(date +%Y-%m-%d)_nextstrain; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

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

###############################################################################
# Merge the metadata with GISAID metadata.-------------------------------------
###############################################################################

# RUN R-SCRIPT. OUTPUTS A FILE "include.txt" for which seqs to include.

###############################################################################
# Merge sequences and QC.------------------------------------------------------
###############################################################################

GISAID_FASTA=$(findTheLatest "${DISTDIR}/global_data/*fasta")

cat $SEQS $GISAID_FASTA | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - > $OUTDIR/seqs.fasta # cat and make oneliner fasta.

# Subset to ones in metadata.
#awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' $OUTDIR/seqs.fasta | # tabularize.
#awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if (">"$1 in seqs) {print seqs[">"$1]} else {print $0" had no matching sequence, excluding." > "/dev/stderr"}}' - $OUTDIR/include.txt | # Subset to ones in metadata.
#tr "\t" "\n" > tmp && mv tmp $OUTDIR/seqs.fasta # de-tabularize.

# Filter bad sequences.
MAXN=5000
MINLENGT=25000

touch $OUTDIR/filtered.txt

cat $OUTDIR/seqs.fasta | 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
awk -v THR=$MAXN -v LEN=$MINLENGTH -v outdir=$OUTDIR '!/^>/ { next } { getline seq; seq2=seq; Nn=gsub(/N/,"",seq) }; {if (length(seq2) > LEN && Nn <= THR) { print $0 "\n" seq2 } else {sub(/^>/,"",$0); print $0 >> outdir"/filtered.txt"}}' - > tmp && mv tmp $OUTDIR/seqs.fasta # Tidy header.

echo "$(wc -l $OUTDIR/filtered.txt | sed 's/ .*//') genomes failed QC and will also be removed from meta.tsv, see $OUTDIR/filtered.txt"

# Subset metadata to remove those sequences. 
#grep -vf $OUTDIR/filtered.txt $OUTDIR/meta.tsv > tmp && mv tmp meta.tsv

# Align and mask bases.
SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
singularity --silent exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "$THISDIR/AlignMask.sh $OUTDIR seqs.fasta $REF/MN908947.3.gb $DEPEND_DIR $THREADS" 

exit 1


###############################################################################
# Run nextstrain - DK only.----------------------------------------------------
###############################################################################

# [RUN FOR ONLY DK SAMPLES AND GET THE CLADE ASSIGNMENT, DUMP "clades.tsv" IN OUTPUT FOLDER]

###############################################################################
# Run nextstrain - global.-----------------------------------------------------
###############################################################################

# [RUN FOR DK + GLOBAL SUBSET, USE "clades.tsv" TO ASSIGN CLADES]

###############################################################################
# Add clades to metadata.------------------------------------------------------
###############################################################################

# [RUN SOME SCRIPT THAT ADDS CLADES TO ALL SEQUENCES IN "seqs.fasta" AND APPENDS TO METADATA]


