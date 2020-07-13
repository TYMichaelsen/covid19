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

#META=/srv/rbd/covid19/genomes/2020-07-08-20-29_export/metadata.tsv
#SEQS=/srv/rbd/covid19/genomes/2020-07-08-20-29_export/sequences.fasta

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

GISAID_META=$(findTheLatest "${DISTDIR}/global_data/*tsv")

Rscript ${THISDIR}/merge_clean_metadata.R -l $META -g $GISAID_META -o $OUTDIR

###############################################################################
# Merge sequences and QC.------------------------------------------------------
###############################################################################

MAXN=5000
MINLENGT=25000

# GISAID data -----------------------------------------------------------------  
GISAID_FASTA=$(findTheLatest "${DISTDIR}/global_data/*fasta")
  
if [ ! -z "$(grep "AlignMasked" <<< $GISAID_FASTA)" ]; then
  echo ""
  echo "GISAID sequences has already been QC'ed, aligned and masked. Skipping that..."
  echo "" 
  
  GISAID_FASTA_OUT=$GISAID_FASTA
  GISAID_FILTER_OUT=$(sed 's/_AlignMasked.fasta/_filtered.txt/' <<< $GISAID_FASTA_OUT)
  
else
  echo ""
  echo "QC'ing, align and mask the GISAID sequences..."
  echo ""
  
  GISAID_FASTA_OUT=$(sed 's/.fasta/_AlignMasked.fasta/' <<< $GISAID_FASTA)
  GISAID_FILTER_OUT=$(sed 's/.fasta/_filtered.txt/' <<< $GISAID_FASTA)
  
  touch $GISAID_FILTER_OUT
  
  # Filter bad sequences.
  cat $GISAID_FASTA | 
  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
  awk -v THR=$MAXN -v LEN=$MINLENGTH -v filtfile=$GISAID_FILTER_OUT '!/^>/ { next } { getline seq; seq2=seq; Nn=gsub(/N/,"",seq) }; {if (length(seq2) > LEN && Nn <= THR) { print $0 "\n" seq2 } else {sub(/^>/,"",$0); print $0 >> filtfile}}' - > $GISAID_FASTA_OUT # Tidy header.
    
  # Align and mask bases.
  SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
  singularity --silent exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "$THISDIR/AlignMask.sh ${DISTDIR}/global_data $GISAID_FASTA_OUT $REF/MN908947.3.gb $DEPEND_DIR $THREADS" 
  
  # Cleanup and rename.
  rm ${DISTDIR}/global_data/aligned*
  rm ${DISTDIR}/global_data/log.out
  mv ${DISTDIR}/global_data/masked.fasta $GISAID_FASTA_OUT
  
fi

echo "$(wc -l $GISAID_FILTER_OUT | sed 's/ .*//') genomes in GISAID data failed QC and will also be removed from metadata, see $OUTDIR/filtered.txt"

# In-house data ---------------------------------------------------------------
echo ""
echo "QC'ing, align and mask the in-house sequences..."
echo ""

# Filter bad sequences.
touch $OUTDIR/filtered.txt

cat $SEQS | 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
awk -v THR=$MAXN -v LEN=$MINLENGTH -v outdir=$OUTDIR '!/^>/ { next } { getline seq; seq2=seq; Nn=gsub(/N/,"",seq) }; {if (length(seq2) > LEN && Nn <= THR) { print $0 "\n" seq2 } else {sub(/^>/,"",$0); print $0 >> outdir"/filtered.txt"}}' - > tmp && mv tmp $OUTDIR/seqs.fasta # Tidy header.

# Align and mask bases.
SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
singularity --silent exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "$THISDIR/AlignMask.sh $OUTDIR $OUTDIR/seqs.fasta $REF/MN908947.3.gb $DEPEND_DIR $THREADS" 

echo "$(wc -l $OUTDIR/filtered.txt | sed 's/ .*//') genomes in in-house data failed QC and will also be removed from metadata, see $OUTDIR/filtered.txt"

# Merge the SSI and GISAID filtered sequences.
cat $GISAID_FILTER_OUT $OUTDIR/filtered.txt > tmp && mv tmp $OUTDIR/filtered.txt

# Merge sequences and filter from metadata ------------------------------------

# Merge GISAID and SSI sequences.
cat $OUTDIR/masked.fasta $GISAID_FASTA_OUT > tmp && mv tmp $OUTDIR/masked.fasta

# Subset to ones in metadata.
awk 'NR > 1 {print $1}' $OUTDIR/metadata_nextstrain.tsv > $OUTDIR/include.txt

cat $OUTDIR/masked.fasta | 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' - | # tabularize.
awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if (">"$1 in seqs) {print seqs[">"$1]} else {print $0" had no matching sequence, excluding." > "/dev/stderr"}}' - $OUTDIR/include.txt | # Subset to ones in metadata.
tr "\t" "\n" > tmp && mv tmp $OUTDIR/masked.fasta # de-tabularize.

# Subset metadata to remove filtered sequences. 
grep -vf $OUTDIR/filtered.txt $OUTDIR/metadata_nextstrain.tsv > tmp && mv tmp $OUTDIR/metadata_nextstrain.tsv

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


