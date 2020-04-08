#!/bin/bash

VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-g file -b file -i string -s] 
-- COVID-19 data generation pipeline v. $VERSION: quantify read mappings in both directions. 

Arguments:
    -h  Show this help text.
    -d  Directory with fastq files.
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
    d) FASTQDIR=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${FASTQDIR+x} ]; then echo "-d $MISSING"; exit 1; fi;
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
REF=/srv/rbd/covid19/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta

mkdir -p $OUTDIR/data
mkdir -p $OUTDIR/tempo/TMPDIR/
mkdir -p $OUTDIR/tempo/trim
#mkdir -p $OUTDIR/tempo/demultiplexed/
mkdir -p $OUTDIR/tempo/articminion/
mkdir -p $OUTDIR/results/genomes/
mkdir -p $OUTDIR/results/coverage/
mkdir -p $OUTDIR/results/mapped_fastq/
mkdir -p $OUTDIR/results/N_counts/

#################################################
# Basecall reads	 							#
#################################################

#################################################
# Pre-process samples 							#
#################################################

#cp $FASTQDIR/*.fastq $OUTDIR/tempo/demultiplexed/

#################################################
# Process samples 								#
#################################################
# Run reference "assemblies"

FILES=$FASTQDIR/*.fastq

for f in $FILES; do

  SAMPLE=$(basename $f | sed 's/.fastq//')
  OUTPUTFILE=$OUTDIR/results/coverage/$SAMPLE.cov.tsv
  
  if [ -s $OUTPUTFILE ]; then 
    echo "$OUTPUTFILE has already been generated"; 
  else

    # Basic artic minion workflow.
    artic minion --medaka --minimap2 --normalise 200 --threads $THREADS --scheme-directory $SCHEMEDIR --read-file $FASTQDIR/$SAMPLE.fastq nCoV-2019/V3.1 $SAMPLE
    mv $SAMPLE* $OUTDIR/tempo/articminion/
    
    # Generate coverage.
    echo -e 'scaffold\tposition\tcoverage' > $OUTDIR/results/coverage/$SAMPLE.cov.tsv
    samtools depth -a $OUTDIR/tempo/articminion/$SAMPLE.sorted.bam >> $OUTDIR/results/coverage/$SAMPLE.cov.tsv
    
    mv $OUTDIR/tempo/articminion/$SAMPLE.consensus.fasta $OUTDIR/results/genomes/
    
    # Extract mapped reads.
    samtools fastq --threads $THREADS -F 4 $OUTDIR/tempo/articminion/$SAMPLE.sorted.bam > $OUTDIR/results/mapped_fastq/$SAMPLE"_virus".fastq
    md5sum $OUTDIR/results/mapped_fastq/$SAMPLE"_virus".fastq >> $OUTDIR/results/md5sums.txt
    
    # Calculate the number of Ns in the sequences.
    awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $OUTDIR/results/genomes/$SAMPLE.consensus.fasta |\
    awk '!/^>/ { next } { getline seq; len=length(seq); Nn=gsub(/N/,"",seq) } {sub(/^>/,"",$0); print $0 "\t" Nn}' - > $OUTDIR/results/N_counts/$SAMPLE"N_count.tsv"
    
  fi
done

#################################################
# Plotting		 								#
#################################################
