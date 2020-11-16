#!/usr/bin/env bash

source activate artic-ncov2019

touch $INDIR/human-filtered.txt

# for each fastq, run sanitizeme.
parallel -j $THREADS '

minimap2 -ax map-ont -t 1 {2} {1} | samtools view --threads 1 -u -f 4 - | samtools fastq --threads - > tmp && mv tmp {1}
echo {1} >> $INDIR/human-filtered.txt

' ::: $INDIR/*fastq ::: $HUMREF
