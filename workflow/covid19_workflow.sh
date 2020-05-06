#!/usr/bin/env bash

# Terminal input --------------------------------------------------------------
INPUT_DIR=$1
SCHEME=$2
RUN_ID=${3:-${INPUT_DIR##*/}}
OUT_DIR=${4:-$INPUT_DIR}
SAMPLE_SHEET=${5:-$INPUT_DIR/sample_sheet.csv}
THREADS=${6:-100}

# Define and format vars ------------------------------------------------------------------
WORKFLOWS_PATH="$(dirname "$(readlink -f "$0")")"
COVID19_PATH=${WORKFLOWS_PATH%/*}
SINGIMG="${COVID19_PATH}/singularity/covid19*.sif"
OUT_DIR=$(readlink -f $OUT_DIR)

echo "[$(date +"%T")] Preparing input and output directories"
if [ -d "${INPUT_DIR}/fastq" ]; then
  echo "Processing fastq files from ${INPUT_DIR}/fastq directory..."
  FASTQ_DIR="${INPUT_DIR}/fastq"
elif [ -d "$INPUT_DIR" ]; then
  echo "Processing fastq files from ${INPUT_DIR} directory..."
  FASTQ_DIR="${INPUT_DIR}"
else
   echo ""
   echo "$INPUT_DIR or ${INPUT_DIR}/fastq does not exist"
   echo "Aborting covid19_workflow..."
   echo ""
   return
fi

if [ ! -d "${OUT_DIR}/sample_sheet.csv" ]; then
  cp $SAMPLE_SHEET ${OUT_DIR}/
fi


# This is the full workflow script.
# Make final_output to dump important stuff.
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/final_output

###############################################################################
# Run demultiplexing.
###############################################################################

if [ -d $OUT_DIR/demultiplexed ]; then 
  echo "demultiplexed data exists, skipping this part."
else 
  singularity \
    exec \
    -B $WORKFLOWS_PATH:$WORKFLOWS_PATH \
    -B $FASTQ_DIR:$FASTQ_DIR \
    -B $OUT_DIR:$OUT_DIR \
    -B $INPUT_DIR:$INPUT_DIR \
    $SINGIMG \
    bash -c "
      $WORKFLOWS_PATH/demultiplex.sh \
        $FASTQ_DIR \
        $RUN_ID \
        $SAMPLE_SHEET \
        $OUT_DIR/demultiplexed \
        $THREADS \
        $SCHEME
      "
fi

###############################################################################
# Run processing
###############################################################################

singularity \
  exec \
  -B $WORKFLOWS_PATH:$WORKFLOWS_PATH \
  -B $INPUT_DIR:$INPUT_DIR \
  -B $OUT_DIR:$OUT_DIR \
  $SINGIMG \
  bash -c "
    source activate artic-ncov2019-medaka;
    # Medaka cannot control mem, need to scale down.
    THREADS_MEDAKA=$((($THREADS+1)/3));
    # Rerun artic only if not existing.
    if [ -d $OUT_DIR/processing/articminion ]; then
      Flag='-a'
    fi;
    $WORKFLOWS_PATH/processing.sh \
      -d $OUT_DIR \
      -s nCoV-2019/$SCHEME \
      -o $OUT_DIR/processing \
      -t \$THREADS_MEDAKA \
      \$Flag;
    retn_code=$?;
    if [ \$retn_code == 1 ]; then echo 'ERROR in processing.sh, exitting.'; exit; fi
    "

###############################################################################
# Run QC
###############################################################################

singularity \
  exec \
  -B $WORKFLOWS_PATH:$WORKFLOWS_PATH \
  -B $INPUT_DIR:$INPUT_DIR \
  -B $OUT_DIR:$OUT_DIR \
  --no-home \
  $SINGIMG \
  bash -c "
    source activate nextstrain
    $WORKFLOWS_PATH/QC.sh \
      -i $OUT_DIR \
      -b $RUN_ID \
      -r $WORKFLOWS_PATH/QC.rmd \
      -t $THREADS
    "

################################################################################
# Sweep important data and put in "output"
################################################################################

cp $OUT_DIR/QC/${RUN_ID}.html $OUT_DIR/final_output/${RUN_ID}.html
cp $OUT_DIR/processing/results/consensus.fasta $OUT_DIR/final_output/consensus.fasta
