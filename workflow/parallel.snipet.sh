#!/usr/bin/env bash
source activate artic-ncov2019-medaka

# for each mapping, run sanitizeme.
cat upload/$DT/tmp_toparallel.txt | parallel -j $THREADS --colsep ':' --bar \
                                             '
# Get the corresponding .bam file.
BAMFILE=$(grep {2}_ {3}/tmp_bamfiles)

if [ -z "$BAMFILE" ]; then
  >&2 echo "warning: .bam file for {1} (LIB-ID: {2}) was not in the artic output." 
else
  # Extract mapped reads and remove human reads with the CDC protocol: https://github.com/CDCgov/SanitizeMe
  samtools fastq --threads 1 -F 4 $BAMFILE 2> /dev/null |\
  minimap2 -ax map-ont {4} - -t 1 2> /dev/null |\
  samtools view --threads 1 -u -f 4 - 2> /dev/null |\
  samtools bam2fq --threads 1 - 2> /dev/null |\
  gzip -c - > {3}/mapped_fastq/"$(sed s/\\//_/g <<< {1})".fastq.gz 2> /dev/null
fi' 

#------------------------------------------------------------------------------
