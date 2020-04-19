#!/bin/bash

VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-g file -b file -i string -s] 
-- COVID-19 data generation pipeline v. $VERSION: quantify read mappings in both directions. 

Arguments:
    -h  Show this help text.
    -d  Library pool directory.
    -o  Output directory.
    -t  Number of threads.

Output:
    To come.

Note: 
This is beta software!
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hd:o:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) POOLDIR=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${POOLDIR+x} ]; then echo "-d $MISSING"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then echo "-d $MISSING"; exit 1; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------
# setup output folders.
if [ "$OUTDIR" != "$PWD" ]; then rm -rf $OUTDIR; mkdir $OUTDIR; fi

# Capture stdout to log file.
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>$OUTDIR/log.out 2>&1

# Settings
#RUNID=CJ023
SCHEMEDIR=/srv/rbd/covid19/software/artic-ncov2019/primer_schemes/
REF=$SCHEMEDIR/nCoV-2019/V3.1/*reference.fasta

mkdir -p $OUTDIR/data
mkdir -p $OUTDIR/tempo/TMPDIR/
mkdir -p $OUTDIR/tempo/trim
mkdir -p $OUTDIR/tempo/articminion/
mkdir -p $OUTDIR/results/genomes/
mkdir -p $OUTDIR/results/coverage/
mkdir -p $OUTDIR/results/mapped_fastq/
mkdir -p $OUTDIR/results/N_counts/

# add metadata to output folder. 
cp $POOLDIR/*sequencing.csv $OUTDIR/

FILES=$POOLDIR/demultiplexed/*.fastq

### Basic artic workflow in parallel.##########################################
parallel -j $THREADS \
'
SAMPLE=$(basename {1} | sed 's/.fastq//');

cd {2}/tempo/articminion;

artic minion \
    --medaka \
    --minimap2 \
    --normalise 200 \
    --threads 1 \
    --scheme-directory {3} \
    --read-file {1} \
    nCoV-2019/V3.1 $SAMPLE;

' ::: $FILES ::: $OUTDIR ::: $SCHEMEDIR

### Misc stuff in parallel. ###################################################
parallel -j $THREADS \
'
SAMPLE=$(basename {1} | sed 's/.fastq//');

# Generate coverage.
echo -e 'scaffold\tposition\tcoverage' > {2}/results/coverage/$SAMPLE.cov.tsv
samtools depth -a {2}/tempo/articminion/$SAMPLE.sorted.bam >> {2}/results/coverage/$SAMPLE.cov.tsv

# Extract mapped reads.
samtools fastq --threads 1 -F 4 {2}/tempo/articminion/$SAMPLE.sorted.bam > {2}/results/mapped_fastq/$SAMPLE"_virus".fastq
md5sum {2}/results/mapped_fastq/$SAMPLE"_virus".fastq >> {2}/results/md5sums.txt

# Calculate the number of Ns in the sequences.
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $OUTDIR/results/genomes/$SAMPLE.consensus.fasta |\
awk '!/^>/ { next } { getline seq; len=length(seq); Nn=gsub(/N/,"",seq) } {sub(/^>/,"",$0); print $0 "\t" Nn}' - > $OUTDIR/results/N_counts/$SAMPLE"N_count.tsv"

' ::: $FILES ::: $OUTDIR
