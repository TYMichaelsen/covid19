#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.2.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-d dir -s dir -o dir -t string -a flag] 
-- COVID-19 data generation pipeline v. $VERSION: quantify read mappings in both directions. 

Arguments:
    -h  Show this help text.
    -d  Library pool directory. MUST BE RELATIVE TO CURRENT WORKING DIRECTORY.
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
mkdir -p $OUTDIR
AAU_COVID19_PATH="$(dirname "$(readlink -f "$0")")"
# Logging
LOG_NAME="$OUTDIR/processing_log_$(date +"%Y-%m-%d_%H-%M").txt"
echo "processing log" >> $LOG_NAME
#echo "AAU COVID-19 revision - $(git -C $AAU_COVID19_PATH rev-parse --short HEAD)" >> $LOG_NAME ##Not working in singularity
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

# Dependencies.
SCHEMEDIR=$AAU_COVID19_PATH/dependencies/primer_schemes                
REF=$AAU_COVID19_PATH/dependencies/ref/MN908947.3.fasta
HUMANREF=$AAU_COVID19_PATH/dependencies/ref/human_g1k_v37.fasta

# setup output folders.
rm -rf $OUTDIR/TMPDIR; mkdir $OUTDIR/TMPDIR/
mkdir -p $OUTDIR/results/
mkdir -p $OUTDIR/results/mapped_fastq

# Setup log.
echo "# input settings (created on $(date))" > $OUTDIR/log.out
echo "-d: "$POOLDIR >> $OUTDIR/log.out
echo "-s: "$SCHEMEVERS >> $OUTDIR/log.out
echo "-o: "$OUTDIR >> $OUTDIR/log.out
echo "-t: "$THREADS >> $OUTDIR/log.out
echo "-a: "$RERUN >> $OUTDIR/log.out

FILES=$POOLDIR/demultiplexed/*.fastq

# For artic, need to make sure we are pointing correctly at files.
nrep=$(echo $OUTDIR/articminion | awk -F'/' '{for (i = 1; i <= NF; i++) a=a"../"; print a}')
FILES_ARTIC=$(for i in $FILES; do echo $nrep$i; done)

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
  
  ' ::: $FILES_ARTIC ::: $OUTDIR ::: $SCHEMEDIR ::: $SCHEMEVERS
   
else

  if [ ! -d "$OUTDIR/articminion" ] || [ ! "$(ls -A $OUTDIR/articminion)" ]; then
    echo "-a was set, but no artic output detected. exitting."
    exit 1
  else
    echo "-a was set, skipping artic workflow."
  fi

fi 

# Capture the libraries that for some reason did not produce consensus.
comm -23 <(ls $FILES | sed 's|.*\/||' | sed 's|.fastq||' | sort) <(ls $OUTDIR/articminion/*.consensus.fasta | sed 's|.*\/||' | sed 's|.consensus.fasta||' | sort) > $OUTDIR/results/failed.txt

echo "$(wc -l $OUTDIR/results/failed.txt | sed 's/ .*//') libraries did not produce consensus, see $OUTDIR/results/failed.txt"

### Misc stuff. ###############################################################
echo "Creating $OUTDIR/results/amplicon_count.tsv"

### Amplicon counts.
echo -e 'library_id\tprimer_id\tstart\tend\tcount' > $OUTDIR/results/amplicon_count.tsv

for i in $FILES; do
  SAMPLE="$(basename $i .fastq)";
  if [ -f $OUTDIR/articminion/$SAMPLE.alignreport.txt ]; then
    awk -v id=$SAMPLE 'NR > 1 { a[$4"\t"$11"\t"$12]++ } END {for (i in a) { print id"\t"i"\t"a[i] }}' $OUTDIR/articminion/$SAMPLE.alignreport.txt >> $OUTDIR/results/amplicon_count.tsv
  fi
done

### Coverage. #################################################################

#echo -e 'library_id\tposition\tcoverage' > $OUTDIR/results/coverage.tsv

# Function for moving average.
#mvavg () {
#  FILE=$1
#  MVINT=$2
#  NAMECOL=$3
#  
#  awk -v id=$NAMECOL -v mvint=$MVINT '{sum+=$3} (NR%10)==0 {print id"\t"$2-(mvint/2)"\t"sum/mvint; sum=0}' $FILE
#}
#export -f mvavg

# Compute coverage.
#parallel -j $THREADS \
#'
#SAMPLE="$(basename {1} .fastq)";

#if [ -f {2}/articminion/$SAMPLE.sorted.bam ]; then
#  samtools depth -a -d 0 {2}/articminion/$SAMPLE.sorted.bam | mvavg - 10 $SAMPLE > {2}/TMPDIR/$SAMPLE.cov.tsv;
#fi
#' ::: $FILES ::: $OUTDIR ::: $HUMANREF

#cat $OUTDIR/TMPDIR/*.cov.tsv >> $OUTDIR/results/coverage.tsv
#rm $OUTDIR/TMPDIR/*.cov.tsv

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
echo -e 'library_id\tposition\tfrac_ALT' > $OUTDIR/results/naive_vcf.tsv

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

echo "Creating $OUTDIR/results/N_counts.tsv"

echo -e 'library_id\tn_count' > $OUTDIR/results/N_counts.tsv;
for i in $FILES; do
  SAMPLE="$(basename $i .fastq)";
  if [ -f $OUTDIR/articminion/$SAMPLE.consensus.fasta ]; then
    awk -v ID=$SAMPLE '!/^>/ { Ns+=gsub(/N/,"")} END {print ID"\t"Ns }' $OUTDIR/articminion/$SAMPLE.consensus.fasta >> $OUTDIR/results/N_counts.tsv;
  fi
done

### Fetch longshot .vcf files #################################################

echo "Creating $OUTDIR/results/longshot_all.tsv"

echo -e "library_id\tposition\tref\talt\tstring" > $OUTDIR/results/longshot_all.tsv
grep -v "#" $OUTDIR/articminion/*longshot.vcf |
awk -F'\t' '{sub(/.longshot.*/,"",$1); sub(/.*\//,"",$1); print $1"\t"$2"\t"$4"\t"$5"\t"$8}' >> $OUTDIR/results/longshot_all.tsv 

### Fetch masked regions ######################################################

echo "Creating $OUTDIR/results/cov_mask_all.tsv"

echo -e "library_id\tstart\tend" > $OUTDIR/results/cov_mask_all.tsv
awk -F'\t' '{sub(/.coverage_mask.txt.*/,"",FILENAME); sub(/.*\//,"",FILENAME); print FILENAME"\t"$2"\t"$3}' $OUTDIR/articminion/*mask.txt >> $OUTDIR/results/cov_mask_all.tsv

### Fetch pass/fail .vcf files ###############################################

echo "Creating $OUTDIR/results/artic_vcf.tsv"

# Positions of SNVs and Ns.
PASSFAIL=$OUTDIR/articminion/*fail.vcf

paste <(echo -e "library_id\ttype") <(grep "^#CHROM" $(echo $PASSFAIL | sed 's/ .*//') | sed 's/#//') --delimiter '\t' > $OUTDIR/results/artic_vcf.tsv

# Specify IN_FAIL directly, for testing.
#IN_FAIL=processing/CJ024/articminion/COV003_LIB-CJ024-75-A-1_UDP0215.fail.vcf	

for IN_FAIL in $PASSFAIL; do	
  
  IN_PASS="$(sed s/fail/pass/ <<< $IN_FAIL)"	
  LIBID="$(basename $IN_FAIL .fail.vcf)"
  
  grep -v "^#" <(gunzip -c $IN_PASS.gz) > $OUTDIR/TMPDIR/pass	  
  grep -v "^#" $IN_FAIL > $OUTDIR/TMPDIR/fail	

awk -F'\t' -v id=$LIBID 'FNR == NR {print id"\tfail\t"$0} FNR != NR {print id"\tpass\t"$0}' $OUTDIR/TMPDIR/fail $OUTDIR/TMPDIR/pass | 	  # If fail print ALT as an N, if pass print the actual ALT.
sort -n -k 2 - >> $OUTDIR/results/artic_vcf.tsv

done
rm $OUTDIR/TMPDIR/pass 	 
rm $OUTDIR/TMPDIR/fail


### Dump genomes passing crude QC filter ######################################
rm -f $OUTDIR/results/filtered.txt
MAXN=3000
MINLENGT=25000

for file in $OUTDIR/articminion/*.consensus.fasta; do echo $file; done
if fgrep ">MN908947.3" $file > /dev/null; then

cat $OUTDIR/articminion/*.consensus.fasta | 
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' - | awk 'NR > 1' - | # make one-line fasta.
awk '/^>/ {sub(/\/ARTIC.*$/,"",$0)}1' - |

awk -v THR=$MAXN -v LEN=$MINLENGTH -v outdir=$OUTDIR '!/^>/ { next } { getline seq; seq2=seq; Nn=gsub(/N/,"",seq) }; {if (length(seq2) > LEN && Nn <= THR) { print $0 "\n" seq2 } else {sub(/^>/,"",$0); print $0 >> outdir"/results/filtered.txt"}}' - > $OUTDIR/results/consensus.fasta # Tidy header.

echo "$(wc -l $OUTDIR/results/filtered.txt | sed 's/ .*//') genomes failed QC, see $OUTDIR/results/filtered.txt"

