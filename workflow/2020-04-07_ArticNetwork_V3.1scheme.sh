#!/bin/bash

VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d dir -o dir -t string -a flag] 
-- COVID-19 data generation pipeline v. $VERSION: quantify read mappings in both directions. 

Arguments:
    -h  Show this help text.
    -d  Library pool directory.
    -o  Output directory.
    -t  Number of threads.
    -a  Flag: set if you don't want to rerun artic.

Output:
    To come.

Note: 
If the output folder already exists, rerunning the code will override everything! Setting the -a flag will preserve the artic output, which is computationally most demanding. 
This is beta software!
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hd:o:t:a' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) POOLDIR=$OPTARG;;
    o) OUTDIR=$OPTARG;;
    t) THREADS=$OPTARG;;
    a) RERUN="true";;
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
if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi

# Capture stdout to log file.
#exec 3>&1 4>&2
#trap 'exec 2>&4 1>&3' 0 1 2 3
#exec 1>$OUTDIR/log.out 2>&1

# Settings
SCHEMEDIR=/srv/rbd/covid19/software/artic-ncov2019/primer_schemes
REF=$SCHEMEDIR/nCoV-2019/V3.1/*reference.fasta

rm -rf $OUTDIR/TMPDIR; mkdir $OUTDIR/TMPDIR/
mkdir -p $OUTDIR/results/
mkdir -p $OUTDIR/results/mapped_fastq/

# add metadata to output folder. 
cp $POOLDIR/*sequencing.csv $OUTDIR/

FILES=$POOLDIR/demultiplexed/*.fastq

### Basic artic workflow in parallel.##########################################
if [ -z ${RERUN+x} ]; then
  mkdir -p $OUTDIR/articminion/

  parallel -j $THREADS \
  '
  SAMPLE=$(basename {1} | sed 's/.fastq//');
  
  cd {2}/articminion;
  
  artic minion \
      --medaka \
      --minimap2 \
      --normalise 200 \
      --threads 1 \
      --scheme-directory {3} \
      --read-file {1} \
      nCoV-2019/V3.1 $SAMPLE;
  
  ' ::: $FILES ::: $OUTDIR ::: $SCHEMEDIR
   
else

  if [ ! -d "$OUTDIR/articminion" ] || [ ! "$(ls -A $OUTDIR/articminion)" ]; then
    echo "-a was set, but no artic output detected. exitting."
    exit 1
  else
    echo "-a was set, skipping artic workflow."
  fi

fi 

### Misc stuff in parallel. ###################################################
echo -e 'LIB_ID\tposition\tcoverage' > $OUTDIR/results/coverage.tsv
touch $OUTDIR/results/md5sums.txt

# Function for moving average.
mvavg () {
  FILE=$1
  MVINT=$2
  NAMECOL=$3
  
  awk -v id=$NAMECOL -v mvint=$MVINT '{sum+=$3} (NR%10)==0 {print id"\t"$2-(mvint/2)"\t"sum/mvint; sum=0}' $FILE
}
export -f mvavg

parallel -j $THREADS \
'
SAMPLE="$(basename {1} .fastq)";

if [ -f {2}/articminion/$SAMPLE.sorted.bam ]; then
  # Compute coverage.
  samtools depth -a -d 0 {2}/articminion/$SAMPLE.sorted.bam | mvavg - 10 $SAMPLE > {2}/TMPDIR/$SAMPLE.cov.tsv;
  # Extract mapped reads.
  samtools fastq --threads 1 -F 4 {2}/articminion/$SAMPLE.sorted.bam > {2}/results/mapped_fastq/$SAMPLE"_virus".fastq;
  md5sum {2}/results/mapped_fastq/$SAMPLE"_virus".fastq > {2}/TMPDIR/$SAMPLE.md5sums.txt;
fi
' ::: $FILES ::: $OUTDIR

cat $OUTDIR/TMPDIR/*.cov.tsv >> $OUTDIR/results/coverage.tsv
rm $OUTDIR/TMPDIR/*.cov.tsv

cat $OUTDIR/TMPDIR/*.md5sums.txt >> $OUTDIR/results/md5sums.txt
rm $OUTDIR/TMPDIR/*.md5sums.txt

### Call naive variants in parallel. ##########################################
MINALT=0.05 # the minimum fraction of bases which supports the alternative base.
BAMFILES=$OUTDIR/articminion/*_1.sorted.bam

# Define function to filter vcf.
filter_vcf () {
  FILE=$1
  CUTOFF=$2
  
  awk -v cutoff=$CUTOFF '/##bcftoolsCommand/ {split($0,a,/ /); sub(/.*\//,"",a[7]); sub(/\.primer.*/,"",a[7]); id=a[7]} !/#/ {
  DP4=$8
  
  # Split into 1) forward REF, 2) reverse REF, 3) forward ALT, 4) reverse ALT. 
  sub(/.*DP4=/,"",DP4)
  sub(/;.*/,"",DP4)
  split(DP4,a,",")
  
  # Compute fraction of ALT.
  if (a[1]+a[2]+a[3]+a[4] == 0)
    fracALT=0
  else
    fracALT=(a[3]+a[4])/(a[1]+a[2]+a[3]+a[4])
  
  # Write .tsv if above filter.
  if (fracALT >= cutoff) print id"\t"$2"\t"fracALT}' $FILE
}
export -f filter_vcf

# Run pileup and filter in parallel.
echo -e 'LIB_ID\tposition\tfrac_ALT' > $OUTDIR/results/naive_vcf.tsv

parallel -j $THREADS \
' 
IN_NAME="$(sed s/\\_1\\.sorted\\.bam// <<< {1})";
OUT_NAME="$(basename $IN_NAME .primertrimmed.nCoV-2019)";

# get base frequency at all variants and prefilter.
bcftools mpileup --skip-indels --threads 1 --fasta-ref {2} $IN_NAME*bam |
bcftools call -c - | filter_vcf - {4} > {3}/TMPDIR/$OUT_NAME.naivevcf.tsv

' ::: $BAMFILES ::: $REF ::: $OUTDIR ::: $MINALT

cat $OUTDIR/TMPDIR/*.naivevcf.tsv >> $OUTDIR/results/naive_vcf.tsv
rm $OUTDIR/TMPDIR/*.naivevcf.tsv

### Count number of Ns ########################################################
echo -e 'LIB_ID\tNs' > $OUTDIR/results/N_counts.tsv;
for i in $FILES; do
  SAMPLE="$(basename $i .fastq)";
  awk -v ID=$SAMPLE '!/^>/ { Ns+=gsub(/N/,"")} END {print ID"\t"Ns }' $OUTDIR/articminion/$SAMPLE.consensus.fasta >> $OUTDIR/results/N_counts.tsv;
done




exit 1


# Get coverage of each amplicon.
i=/srv/rbd/covid19/data/CJ023/LIB-CJ023-POOL1/demultiplexed/COV001_LIB-CJ023-1-A-1_UDP0377.fastq

SAMPLE="$(basename $i .fastq)";
align_trim $SCHEMEDIR/nCoV-2019/V3.1/nCoV-2019.scheme.bed --remove-incorrect-pairs < $OUTDIR/articminion/$SAMPLE.sorted.bam |
samtools sort -T $SAMPLE - -o $SAMPLE.primertrimmed.nonorm.sorted.bam 

samtools index $SAMPLE.primertrimmed.nonorm.sorted.bam 

# Split into pool 1 and 2.
samtools view -b -r "nCoV-2019_2" $SAMPLE.primertrimmed.nonorm.sorted.bam > $SAMPLE.primertrimmed.nonorm.2.sorted.bam
samtools index $SAMPLE.primertrimmed.nonorm.2.sorted.bam
samtools view -b -r "nCoV-2019_1" $SAMPLE.primertrimmed.nonorm.sorted.bam > $SAMPLE.primertrimmed.nonorm.1.sorted.bam
samtools index $SAMPLE.primertrimmed.nonorm.1.sorted.bam

# Calculate average coverage.
