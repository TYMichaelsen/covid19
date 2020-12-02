#!/usr/bin/env bash

source activate artic-ncov2019

mkdir $OUTdir

# Logging
LOG_NAME="$OUTdir/human-filtering_log_$(date +"%Y-%m-%d_%H-%M").txt"
echo "human read filtering log" >> $LOG_NAME
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

# for each fastq, map to human reference and remove reads that mapped.
parallel -j $THREADS '

minimap2 -ax map-ont -t 1 {3} {1} | samtools view --threads 1 -u -f 4 - | samtools fastq --threads - > {2}/{1/}

' ::: $INdir/*fastq ::: $OUTdir ::: $HUMREF