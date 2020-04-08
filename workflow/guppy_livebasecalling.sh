#!/bin/bash 

# Terminal input
RUN_DIR=$1
FASTQ_DIR=${2:-fastq}
THREADS=${3:-10}
CONFIG=${3:-dna_r9.4.1_450bps_hac.cfg}

# Preparation -----------------------------------------------------------------
echo "[$(date +"%T")] Preparing basecalling session.."

# Check if starting or continuing basecalling session
if [ ! -d "$FASTQ_DIR" ]; then
  mkdir "$FASTQ_DIR"
  echo "#Basecalled fast5 files" > $FASTQ_DIR/basecalled_files.txt
else
  echo "$FASTQ_DIR exists... "
  # Check if basecalled_files.txt exists
  if [ ! -f "$FASTQ_DIR/basecalled_files.txt" ]; then
    # No - stop basecalling
    echo "$FASTQ_DIR/basecalled_files.txt does not exist."
    echo "Can't continue basecalling session safely."
    echo "Clean up output folder and restart basecalling."
    echo "Exiting..."
    exit 1
  else
    # Yes - continue basecalling
    echo "$(wc -l < $FASTQ_DIR/basecalled_files.txt) fast5 files processed in session already."
    echo "Continuing basecalling session..."
    BATCH=$(find $FASTQ_DIR -type f |\
      sed -e "s|$FASTQ_DIR/||" -e 's/_.*//'|\
      sort -n |\
      tail -n 1)
    rm -rf $FASTQ_DIR/processing/
  fi
fi

# Basecalling -----------------------------------------------------------------
echo "[$(date +"%T")] Starting live basecalling loop..."
while : ; do
  
  # Find new fast5 files
  FAST5_NEW=$(
    find "$RUN_DIR" -name "*.fast5" -type f |\
    grep -vFf $FASTQ_DIR/basecalled_files.txt
  )

  FAST5_NEW_N=$(echo "$FAST5_NEW" | grep -c '[^[:space:]]')

  # Start basecalling if new fast5 files was found
  if [ "$FAST5_NEW_N" -eq 0 ]; then
    # No new fast5 files was found
    ATTEMPT=$((ATTEMPT+1))
    if [ "$ATTEMPT" -eq 15 ]; then
      # No new fast5 files for 15 min. Stop basecalling
      echo "No new fast5 files detected in $RUN_DIR for 15 minutes"
      echo "Stopping basecalling..."
      exit 1
    fi
  else
    # New fast5 files was found
    echo "[$(date +"%T")] Starting basecalling of $FAST5_NEW_N fast5 files ..."

    # Update counters
    BATCH=$((BATCH+1))
    ATTEMPT=0

    # Basecall fast5 batch
    guppy_basecaller \
      --save_path $FASTQ_DIR/processing \
      -c $CONFIG \
      --input_file_list <(echo "$FAST5_NEW") \
      --device 'auto'
  
    # Post basecalling
    if [ -f "$FASTQ_DIR/processing/sequencing_telemetry.js" ]; then
      # Register basecalled files
      echo "$FAST5_NEW" >> $FASTQ_DIR/basecalled_files.txt

      # Transfer files
      for FILE in $FASTQ_DIR/processing/*; do
        mv $FILE $FASTQ_DIR/${BATCH}_${FILE##*/}
      done
    fi
    rm -rf $FASTQ_DIR/processing/
  fi

  # Wait for new fast5 files
  if [ "$ATTEMPT" -gt 0 ]; then
    echo "[$(date +"%T")] Waited $ATTEMPT minute(s) for new fast5 files. Check again in 1 minute. "
    echo "Press [CTRL+C] to stop basecalling loop."
    sleep 60 #sec
  fi
done

# Testing
#find "$RUN_DIR" -name "*fast5" | head -n 10 | xargs cp -t fast5/