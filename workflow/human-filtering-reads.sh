#!/usr/bin/env bash

source activate artic-ncov2019

# for each fastq, run sanitizeme.
parallel -j $THREADS '

minimap2 -ax map-ont -t 1 {2} {1}.fastq | samtools view --threads 1 -u -f 4 - | samtools fastq --threads - > {1}.humfilt.fastq

' ::: `ls $INDIR/*fastq | sed 's/.fastq//'` ::: $HUMREF
