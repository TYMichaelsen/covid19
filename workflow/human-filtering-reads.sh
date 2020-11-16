#!/usr/bin/env bash

source activate artic-ncov2019

# for each fastq, run sanitizeme.
parallel -j $THREADS '

minimap2 -ax map-ont -t 1 {2} {1} | samtools view --threads 1 -u -f 4 - | samtools fastq --threads - > {1}.tmp && mv {1}.tmp {1}

' ::: $INDIR/*fastq ::: $HUMREF

# Save that the files has been filtered.
cat $INDIR/*fastq.txt > $INDIR/human-filtered.txt