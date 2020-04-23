#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.2.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d dir -s dir -o dir -t string -a flag] 
-- COVID-19 data generation pipeline v. $VERSION: quantify read mappings in both directions. 

Arguments:
    -h  Show this help text.
    -d  Library pool directory.
    -s  Scheme version directory.
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
while getopts ':hd:s:o:t:a' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    d) POOLDIR=$OPTARG;;
    s) SCHEMEVERS=$OPTARG;;
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
if [ -z ${SCHEMEVERS+x} ]; then echo "-s $MISSING"; exit 1; fi;
if [ -z ${OUTDIR+x} ]; then echo "-d $MISSING"; exit 1; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------
# setup output folders.
if [ ! -d "$OUTDIR" ]; then mkdir $OUTDIR; fi

# Capture stdout to log file.
#exec 3>&1 4>&2
#trap 'exec 2>&4 1>&3' 0 1 2 3
#exec 1>$OUTDIR/log.out 2>&1

# Setup log.
echo "# input settings" > $OUTDIR/log.out
echo "-d: "$POOLDIR >> $OUTDIR/log.out
echo "-s: "$SCHEMEVERS >> $OUTDIR/log.out
echo "-o: "$OUTDIR >> $OUTDIR/log.out
echo "-t: "$THREADS >> $OUTDIR/log.out
echo "-a: "$RERUN >> $OUTDIR/log.out

# Settings
SCHEMEDIR=/srv/rbd/covid19/git/covid19/workflow/primer_schemes
REF=/srv/rbd/covid19/data/reference/MN908947.3.fasta

rm -rf $OUTDIR/TMPDIR; mkdir $OUTDIR/TMPDIR/
mkdir -p $OUTDIR/results/

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
      {4} $SAMPLE;
  
  ' ::: $FILES ::: $OUTDIR ::: $SCHEMEDIR ::: $SCHEMEVERS
   
else

  if [ ! -d "$OUTDIR/articminion" ] || [ ! "$(ls -A $OUTDIR/articminion)" ]; then
    echo "-a was set, but no artic output detected. exitting."
    exit 1
  else
    echo "-a was set, skipping artic workflow."
  fi

fi 

### Misc stuff. ###############################################################
echo "Creating $OUTDIR/results/amplicon_count.tsv"

### Amplicon counts.
echo -e 'LIB_ID\tprimer_id\tstart\tend\tcount' > $OUTDIR/results/amplicon_count.tsv

for i in $FILES; do
  SAMPLE="$(basename $i .fastq)";
  if [ -f $OUTDIR/articminion/$SAMPLE.alignreport.txt ]; then
    awk -v id=$SAMPLE 'NR > 1 { a[$4"\t"$11"\t"$12]++ } END {for (i in a) { print id"\t"i"\t"a[i] }}' $OUTDIR/articminion/$SAMPLE.alignreport.txt >> $OUTDIR/results/amplicon_count.tsv
  fi
done

echo -e 'LIB_ID\tposition\tcoverage' > $OUTDIR/results/coverage.tsv
touch $OUTDIR/results/md5sums.txt

### Coverage.
# Function for moving average.
mvavg () {
  FILE=$1
  MVINT=$2
  NAMECOL=$3
  
  awk -v id=$NAMECOL -v mvint=$MVINT '{sum+=$3} (NR%10)==0 {print id"\t"$2-(mvint/2)"\t"sum/mvint; sum=0}' $FILE
}
export -f mvavg

# Compute coverage.
parallel -j $THREADS \
'
SAMPLE="$(basename {1} .fastq)";

if [ -f {2}/articminion/$SAMPLE.sorted.bam ]; then
  samtools depth -a -d 0 {2}/articminion/$SAMPLE.sorted.bam | mvavg - 10 $SAMPLE > {2}/TMPDIR/$SAMPLE.cov.tsv;
fi
' ::: $FILES ::: $OUTDIR ::: $HUMANREF

cat $OUTDIR/TMPDIR/*.cov.tsv >> $OUTDIR/results/coverage.tsv
rm $OUTDIR/TMPDIR/*.cov.tsv

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

echo "Creating $OUTDIR/results/naive_vcf.tsv"

parallel -j $THREADS \
' 
IN_NAME="$(sed s/\\_1\\.sorted\\.bam// <<< {1})";
OUT_NAME="$(basename $IN_NAME .primertrimmed.nCoV-2019)";

# get base frequency at all variants and prefilter.
bcftools mpileup --skip-indels --threads 1 --fasta-ref {2} $IN_NAME*bam 2> /dev/null |
bcftools call -c - 2> /dev/null | filter_vcf - {4} > {3}/TMPDIR/$OUT_NAME.naivevcf.tsv

' ::: $BAMFILES ::: $REF ::: $OUTDIR ::: $MINALT

cat $OUTDIR/TMPDIR/*.naivevcf.tsv >> $OUTDIR/results/naive_vcf.tsv
rm $OUTDIR/TMPDIR/*.naivevcf.tsv

### Count number of Ns ########################################################
echo -e 'LIB_ID\tNs' > $OUTDIR/results/N_counts.tsv;
for i in $FILES; do
  SAMPLE="$(basename $i .fastq)";
  if [ -f $OUTDIR/articminion/$SAMPLE.consensus.fasta ]; then
    awk -v ID=$SAMPLE '!/^>/ { Ns+=gsub(/N/,"")} END {print ID"\t"Ns }' $OUTDIR/articminion/$SAMPLE.consensus.fasta >> $OUTDIR/results/N_counts.tsv;
  fi
done

### Dump genomes passing crude QC filter ######################################
rm -f $OUTDIR/results/filtered.txt
MAXN=5000
MINLENGT=25000

cat $OUTDIR/articminion/*.consensus.fasta | 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
awk '/^>/ {sub(/_UDP.*$/,"",$0);sub(/.*_LIB/,">LIB",$0)}1' - |
awk -v THR=$MAXN -v LEN=$MINLENGTH -v outdir=$OUTDIR '!/^>/ { next } { getline seq; seq2=seq; Nn=gsub(/N/,"",seq) }; {if (length(seq2) > LEN && Nn <= THR) { print $0 "\n" seq2 } else {sub(/^>/,"",$0); print $0 >> outdir"/results/filtered.txt"}}' - > $OUTDIR/results/consensus.fasta # Tidy header.

echo $(wc -l $OUTDIR/results/filtered.txt | sed 's/ .*//') genomes failed QC, see $OUTDIR/results/filtered.txt

exit 1

# Positions of SNVs and Ns.
echo -e 'LIB_ID\tposition\tALT' > $OUTDIR/results/artic_vcf.tsv

PASSFAIL=$OUTDIR/articminion/*fail.vcf

IN_FAIL=CJ024/articminion/COV003_LIB-CJ024-75-A-1_UDP0215.fail.vcf

for IN_FAIL in $PASSFAIL; do
  IN_PASS="$(sed s/fail/pass/ <<< $IN_FAIL)"
  LIBID="$(basename $IN_FAIL .fail.vcf)"
  
  grep -v "^#" <(gunzip -c $IN_PASS.gz) > $OUTDIR/TMPDIR/pass
  grep -v "^#" $IN_FAIL > $OUTDIR/TMPDIR/fail
  
  awk -F'\t' -v id=$LIBID 'FNR == NR {print id"\t"$2"\tN"} FNR != NR {print id"\t"$2"\t"$5}' $OUTDIR/TMPDIR/fail $OUTDIR/TMPDIR/pass | 
  sort -n -k 2 - >> $OUTDIR/results/artic_vcf.tsv
  
  rm $OUTDIR/TMPDIR/pass 
  rm $OUTDIR/TMPDIR/fail
done