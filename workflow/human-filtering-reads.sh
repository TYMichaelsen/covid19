#!/usr/bin/env bash

source activate artic-ncov2019-medaka

# Logging
LOG_NAME="$INDIR/filtering_log_$(date +"%Y-%m-%d_%H-%M").txt"
echo "filtering log" >> $LOG_NAME
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

# for each fastq, run sanitizeme.
parallel -j $THREADS '

minimap2 -ax map-ont -t 1 {2} {1} | samtools view --threads 1 -u -f 4 - | samtools fastq --threads 1 - > {1}.tmp && mv {1}.tmp {1}

' ::: $INDIR/*fastq ::: $HUMREF

# Save that the files has been filtered.
echo "Data has been human filtered" > $INDIR/filtered.txt
