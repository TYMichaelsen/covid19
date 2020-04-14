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
# NB: 20 bp of each terminal adaptor was targeted
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
echo ""
echo "[$(date +"%T")] Demultiplexing trimmed reads"
echo ""
cutadapt \
  -e 0.1 \
  -O 10 \
  -m 450 \
  -M 1800 \
  -g file:$OUT_DIR/barcodes_used.fasta \
  -o "$OUT_DIR/{name}.tmp" \
  $OUT_DIR/reads_trim.fq

# Trim adapter and primers
echo ""
echo "[$(date +"%T")] Trimming remaining adapter sequence"
echo ""
BARCODE_LIST=$(gawk -F "," 'NR > 1{print $18}' $METADATA)
for BARCODE in $BARCODE_LIST; do
  BARCODE_FILE=$(find $OUT_DIR -name "*$BARCODE*tmp")
  BARCODE_BASE=${BARCODE_FILE%.*}
  cutadapt \
    -m 450 \
    -M 1800 \
    -u 15 \
    -u -14 \
    -o ${BARCODE_BASE}.fastq \
    $BARCODE_FILE
done

# Cleanup
echo ""
echo "[$(date +"%T")] Cleaning up temp files"
echo ""
#rm ./reads_trim.fq ./*tmp

echo ""
echo "[$(date +"%T")] Demultiplexing done..."
echo ""
