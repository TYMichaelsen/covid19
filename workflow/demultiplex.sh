#!/bin/bash 

# General settings
FASTQ_DIR=$1
METADATA=$2
OUT_DIR=${3:-demultiplexed}
THREADS=${4:-1}

AAU_COVID19_PATH="$(dirname "$(readlink -f "$0")")"

# Preparation
mkdir $OUT_DIR

# Logging
LOG_NAME="$OUT_DIR/demultiplex_log_$(date +"%Y-%m-%d-%T").txt"
echo "demultiplex log" >> $LOG_NAME
echo "AAU COVID-19 revision - $(git -C $AAU_COVID19_PATH rev-parse --short HEAD)" >> $LOG_NAME
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

# Trim end adaptors and filter by length
#-- 20 bp of each terminal adaptor is targeted
#-- Sequences are reverse complemented to obtain same i7i5 adapter orientation 
echo ""
echo "[$(date +"%T")] Trimming adapter terminals and filtering by length"
echo ""
find \
  $FASTQ_DIR/ \
  -name "*.fastq" \
  -exec cat {} + |\
cutadapt \
  -j $THREADS \
  -e 0.20 \
  -O 10 \
  -m 500 \
  -M 1800 \
  --discard-untrimmed \
  --revcomp \
  -g CAGAAGACGGCATACGAGAT...GTGTAGATCTCGGTGGTCGC \
  -o $OUT_DIR/reads_trim.fq \
  -

# Compile barcode file
echo ""
echo "[$(date +"%T")] Compiling used barcode file from $METADATA"
echo ""
gawk \
  -F "," \
  '
    FNR==NR && FNR > 1 {  
      LIB_NAME[$18]=">" $1 "_" $7 "_" $18
      next
    }
    {
      if ($1 in LIB_NAME){
        print LIB_NAME[$1]
        print $2
      } 
    }
  ' \
  $METADATA \
  $AAU_COVID19_PATH/barcodes.csv \
  > $OUT_DIR/barcodes_used.fasta

# Demultiplex based
#-- Demultiplexing depends on the i7i5 adaptor orientation
echo ""
echo "[$(date +"%T")] Demultiplexing trimmed reads"
echo ""

demux_wrap(){
  # Input
  local BARCODE_FILE=$1
  local OUT_DIR=$2
  
  # Create outdir
  mkdir $OUT_DIR
  
  # Demultiplex
  cutadapt \
    -e 0.2 \
    -O 10 \
    -m 450 \
    -M 1800 \
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
      $OUT_DIR/demux{#}
    "

# Trim internal adapter sequences and revert reverse complement
#-- Revert sequence orientation to prepare for downstream consensus calling
echo ""
echo "[$(date +"%T")] Trimming remaining adapter sequence"
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
  gawk '
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
      getline; SEQ=substr($0, 16, length($0) - 15 - 14)
      getline; SPACER=$0
      getline; QUAL=substr($0, 16, length($0) - 15 - 14)
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
  'NR > 1{print $18}' \
  $METADATA |\
parallel \
  --env trim_revertrc \
  -j $THREADS \
  "trim_revertrc $OUT_DIR {}"

# Cleanup
echo ""
echo "[$(date +"%T")] Cleaning up temp files"
echo ""
#rm -rf $OUT_DIR/reads_trim.fq $OUT_DIR/*tmp

echo ""
echo "[$(date +"%T")] Demultiplexing done..."
echo ""
