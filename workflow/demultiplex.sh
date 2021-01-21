#!/bin/bash 

# General settings
FASTQ_DIR=$1
RUN_ID=$2
METADATA=$3
OUT_DIR=${4:-demultiplexed}
THREADS=${5:-1}
ARTIC_SCHEME=${6:-aau_long_v3.1}

AAU_COVID19_PATH="$(dirname "$(readlink -f "$0")")"

# Preparation
mkdir $OUT_DIR


# Setup ARTIC scheme settings
#-- The workflow expects barcoded adaptors in the following format
#-- <ADP1_outer> <barcode1> <ADP1_inner> <Target sequence> <ADP2_inner> <barcode2> <ADP2_outer>

if [ "$ARTIC_SCHEME" == "v1" ] || [ "$ARTIC_SCHEME" == "v2" ] || [ "$ARTIC_SCHEME" == "v3" ]; then
  AMP_MIN_LENGTH=350
  AMP_MAX_LENGTH=600
  # Expects ONT native barcodes 1-24
  ADP12_OUTER="TTACGTATTGCTAAGGTTAA...TTAACCTTAGCAATACGTAA"
  ADP1_INNER_LEN=8
  ADP2_INNER_LEN=8
elif [ "$ARTIC_SCHEME" == "aau_long_v3.1" ]; then
  AMP_MIN_LENGTH=500
  AMP_MAX_LENGTH=1800
  # Expects ILM nextera type barcodes
  ADP12_OUTER="CAGAAGACGGCATACGAGAT...GTGTAGATCTCGGTGGTCGC"
  ADP1_INNER_LEN=15 
  ADP2_INNER_LEN=15
elif [ "$ARTIC_SCHEME" == "aau_short_v3" ]; then
  AMP_MIN_LENGTH=350
  AMP_MAX_LENGTH=600
  # Expects ILM nextera type barcodes
  ADP12_OUTER="CAGAAGACGGCATACGAGAT...GTGTAGATCTCGGTGGTCGC"
  ADP1_INNER_LEN=15 
  ADP2_INNER_LEN=15
else
  echo "$ARTIC_SCHEME is not a valid ARTIC scheme."
  echo "Available ARTIC schemes are [v1, v2, v3] (Standard ARTIC) or [aau_long_v3.1, aau_short_v3] (AAU ARTIC)."
  echo "Exiting ..."
  exit 1
fi


# Logging
LOG_NAME="$OUT_DIR/demultiplex_log_$(date +"%Y-%m-%d_%H-%M").txt"
echo "demultiplex log" >> $LOG_NAME
# echo "AAU COVID-19 revision - $(git init; git -C $AAU_COVID19_PATH rev-parse --short HEAD)" >> $LOG_NAME ##Not working in singularity
echo "Command: $0 $*" >> $LOG_NAME
exec 1>>$LOG_NAME
exec 2>&1

# Compile barcode file
rm -f $OUT_DIR/missing.txt
touch $OUT_DIR/missing.txt

gawk \
  -F "," \
  -v RUN_ID="$RUN_ID" \
  -v OUT_DIR="$OUT_DIR" \
  '
    FNR==NR {  
      BARCODE[$1]=$2
      next 
    }
    {
      if ($2 in BARCODE){
        print ">" RUN_ID "_" $1 "_" $2
        print BARCODE[$2]
      }   
      else 
      {
        print $2 >> OUT_DIR"/missing.txt"
      }
    }
  ' \
  $AAU_COVID19_PATH/dependencies/demultiplex/barcodes.csv \
  $METADATA \
  > $OUT_DIR/barcodes_used.fasta
 
if [ -s $OUT_DIR/missing.txt ]; then 
  exit 123
else
  echo "All barcodes in $METADATA were legit, continue workflow..."
fi

  
# Trim end adaptors and filter by length
#-- 20 bp of each terminal adaptor is targeted
#-- Outer adaptors are searched for together with cutadapt <ADP1_outer>...<ADP2_outer>
#-- Sequences are reverse complemented to obtain same i7i5 adapter orientation 
echo ""
echo "[$(date +"%T")] Trimming outer adapter sequences and filter by length"
echo ""
find \
  $FASTQ_DIR/ \
  -mindepth 1 \
  -maxdepth 1 \
  -name "*.fastq" \
  -exec cat {} + |\
cutadapt \
  -j $THREADS \
  -e 0.20 \
  -O 10 \
  -m $AMP_MIN_LENGTH \
  -M $AMP_MAX_LENGTH \
  --discard-untrimmed \
  --revcomp \
  -g $ADP12_OUTER \
  -o $OUT_DIR/reads_trim.fq \
  -

# Demultiplex based
#-- Dual barcode demultiplexing requires i7i5 adaptor orientation
#-- If orientation is not streamlined reads might be wrongly assigned as
#-- some dual barcoding schemes re-use barcode sequences in reverse complement orientation
echo ""
echo "[$(date +"%T")] Demultiplexing trimmed reads"
echo ""

demux_wrap(){
  # Input
  local BARCODE_FILE=$1
  local OUT_DIR=$2
  local AMP_MIN_LENGTH=$3
  local AMP_MAX_LENGTH=$4
  
  # Create outdir
  mkdir $OUT_DIR
  
  # Demultiplex
  cutadapt \
    -e 0.2 \
    -O 10 \
    -m $AMP_MIN_LENGTH \
    -M $AMP_MAX_LENGTH \
    -g file:$BARCODE_FILE \
    -o "$OUT_DIR/{name}.tmp" \
    -
}

export -f demux_wrap

cat $OUT_DIR/reads_trim.fq |\
  parallel \
    --env demux_wrap \
    -L4 \
    -j $THREADS \
    --block 300M \
    --pipe \
    "demux_wrap \
      $OUT_DIR/barcodes_used.fasta \
      $OUT_DIR/demux{#} \
      $AMP_MIN_LENGTH \
      $AMP_MAX_LENGTH
    "

# Trim internal adapter sequences and revert reverse complement
#-- Trim of inner adaptors based on adaptor length
#-- Revert sequence orientation to prepare for downstream consensus calling

echo ""
echo "[$(date +"%T")] Trimming inner adapter sequences"
echo ""
trim_revertrc() {

  # Input
  local OUT_DIR=$1
  local BARCODE=$2
  
  # Find file
  BARCODE_FILES=$(find $OUT_DIR -name "*$BARCODE*tmp")
  BARCODE_NAME=$(
    echo $BARCODE_FILES |\
    sed -e 's| .*||g' -e 's|.*/||g' -e 's|.tmp$||'
  )
  
  # Trim adaptor and revert rc
  cat \
    $BARCODE_FILES |\
  gawk \
    -v ADP1_INNER_LEN="$ADP1_INNER_LEN" \
    -v ADP2_INNER_LEN="$ADP2_INNER_LEN" \
    '
    # Define function for reverse complement
    function revcomp(s,  i, o) {
      o = ""
      for(i = length; i > 0; i--)
           o = o c[substr(s, i, 1)]
      return(o)
    }
    function rev(s,  i, o) {
      o = ""
      for(i = length; i > 0; i--)
           o = o substr(s, i, 1)
      return(o)
    }
    BEGIN{
      # Define revcomp mapping vector
      c["A"] = "T"
      c["C"] = "G"
      c["G"] = "C"
      c["T"] = "A"
      
    }
    NR%4==1{
      # Trim record
      HEAD=$0
      getline; SEQ=substr($0, ADP1_INNER_LEN + 1, length($0) - ADP1_INNER_LEN - ADP2_INNER_LEN)
      getline; SPACER=$0
      getline; QUAL=substr($0, ADP1_INNER_LEN + 1, length($0) - ADP1_INNER_LEN - ADP2_INNER_LEN)
      # Print sequence based on orientation
      if(HEAD ~ /.* rc$/){
        sub(/ rc$/, "", HEAD)
        print HEAD "\n" revcomp(SEQ) "\n" SPACER "\n" rev(QUAL)
      } else {
        print HEAD "\n" SEQ "\n" SPACER "\n" QUAL
      }
    }
  ' \
  - \
  > ${OUT_DIR}/${BARCODE_NAME}.fastq
}
export -f trim_revertrc

gawk \
  -F "," \
  '{print $2}' \
  $METADATA |\
parallel \
  --env trim_revertrc \
  -j $THREADS \
  "trim_revertrc $OUT_DIR {}"

# Cleanup
echo ""
echo "[$(date +"%T")] Cleaning up temp files"
echo ""
rm -rf $OUT_DIR/reads_trim.fq $OUT_DIR/demux*

echo ""
echo "[$(date +"%T")] Demultiplexing done..."
echo ""
