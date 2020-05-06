#!/usr/bin/env bash

# Terminal input --------------------------------------------------------------
INPUT_DIR=$1
SCHEME=$2
RUN_ID=${3:-${INPUT_DIR##*/}}
OUT_DIR=${4:-$INPUT_DIR}
SAMPLE_SHEET=${5:-$INPUT_DIR/sample_sheet.csv}
THREADS=${6:-100}

# Preparation -----------------------------------------------------------------

echo ""
echo "[$(date +"%T")] Performing pre-run preparation and checks"
echo ""

# Script paths
WORKFLOW_PATH="$(dirname "$(readlink -f "$0")")"
COVID19_PATH=${WORKFLOW_PATH%/*}
SINGIMG="${COVID19_PATH}/singularity/covid19*.sif"

# Determien absolute path to output
OUT_DIR=$(readlink -f $OUT_DIR)

# Create output folder
mkdir -p $OUT_DIR
mkdir -p $OUT_DIR/final_output

# Copy sample_sheet.csv to output folder
echo $SAMPLE_SHEET
if [ ! -d "${OUT_DIR}/sample_sheet.csv" ]; then
  cp $SAMPLE_SHEET ${OUT_DIR}/
fi

# Check references 
if [ ! -f "${WORKFLOW_PATH}/dependencies/ref/human_g1k_v37.fasta" ]; then
  echo ""
  echo "Downloading human reference genome to ${WORKFLOW_PATH}/dependencies/ref/human_g1k_v37.fasta"
  echo ""
  
  wget \
    -O- \
    'ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz' |\
  gunzip > "${WORKFLOW_PATH}/dependencies/ref/human_g1k_v37.fasta"
fi


# This is the full workflow script.
# Make final_output to dump important stuff.

###############################################################################
# Demultiplexing
###############################################################################

echo ""
echo "[$(date +"%T")] Demultiplexing and trimming of fastq data"
echo ""

if [ -d $OUT_DIR/demultiplexed ]; then 
  echo "$OUT_DIR/demultiplexed directory exists. Demultiplexing is skipped..."
else 
  # Check if default folder structure is present
  if [ -d "${INPUT_DIR}/fastq" ]; then
    echo ""
    echo "Processing fastq files from ${INPUT_DIR}/fastq directory..."
    FASTQ_DIR="${INPUT_DIR}/fastq"
    echo ""
  elif [ -d "$INPUT_DIR" ]; then
    echo ""
    echo "Processing fastq files from ${INPUT_DIR} directory..."
    FASTQ_DIR="${INPUT_DIR}"
    echo ""
  else
    echo ""
    echo "$INPUT_DIR or ${INPUT_DIR}/fastq does not exist"
    echo "Aborting covid19_workflow..."
    echo ""
    return
  fi
  
  # Demultiplexing
  singularity \
    exec \
    -B $WORKFLOW_PATH:$WORKFLOW_PATH \
    -B $FASTQ_DIR:$FASTQ_DIR \
    -B $OUT_DIR:$OUT_DIR \
    -B $INPUT_DIR:$INPUT_DIR \
    $SINGIMG \
    bash -c "
      $WORKFLOW_PATH/demultiplex.sh \
        $FASTQ_DIR \
        $RUN_ID \
        $SAMPLE_SHEET \
        $OUT_DIR/demultiplexed \
        $THREADS \
        $SCHEME
      "
fi


###############################################################################
# Generate genomes
###############################################################################

echo ""
echo "[$(date +"%T")] Generating genomes with ARTIC medaka pipeline"
echo ""

singularity \
  exec \
  -B $WORKFLOW_PATH:$WORKFLOW_PATH \
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
    $WORKFLOW_PATH/processing.sh \
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

echo ""
echo "[$(date +"%T")] Generating genomes with ARTIC medaka pipeline"
echo ""

singularity \
  exec \
  -B $WORKFLOW_PATH:$WORKFLOW_PATH \
  -B $INPUT_DIR:$INPUT_DIR \
  -B $OUT_DIR:$OUT_DIR \
  --no-home \
  $SINGIMG \
  bash -c "
    source activate nextstrain
    $WORKFLOW_PATH/QC.sh \
      -i $OUT_DIR \
      -b $RUN_ID \
      -r $WORKFLOW_PATH/QC.rmd \
      -t $THREADS
    "

################################################################################
# Sweep important data and put in "output"
################################################################################

cp $OUT_DIR/QC/${RUN_ID}.html $OUT_DIR/final_output/${RUN_ID}.html
cp $OUT_DIR/processing/results/consensus.fasta $OUT_DIR/final_output/consensus.fasta
