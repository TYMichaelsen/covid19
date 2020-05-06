#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-s file -m file -t int] 
-- COVID-19 data export for GISAID v. $VERSION:  

Arguments:
    -h  Show this help text.
    -m  Mapping between library_id and user-specified id. Samples not in library_id will be excluded.
    -t  Number of threads.
    -x  Flag: set if you want to also dump fastq reads.

Output:
    To come.

Note: 
    (DEPRECIATED) Using the option -H it is possible to specify an extended sample name that should be used instead of just the ID. 
    Example: 
    ID = HH-1
    -H = "DK/ALAB-@/2020" will lead to renaming as 'DK_ALAB-SSI-HH-1_2020'
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hm:t:x' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) MAP_ID=$OPTARG;;
    t) THREADS=$OPTARG;;
    x) DUMPREADS="true";;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${MAP_ID+x} ]; then echo "-m $MISSING"; exit 1; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

### Code.----------------------------------------------------------------------

GITDIR="$(dirname "$(readlink -f "$0")")"

# How the header should look like.
NAMESTRING="DK/ALAB-@/2020"

# Set dependent directories/data.
#ARTICDIR=( CJ*/processing )
#HUMANREF=$AAU_COVID19_PATH/human_g1k_v37.fasta

# Make time-stamped export subfolder.
DT=$(date +%Y_%m_%d_%H-%M)"_export"
mkdir -p genomes
mkdir -p genomes/$DT

# Sweep batch-folders for consensus sequences (OBS: renaming headers according to AAU library naming scheme).
cat processing/CJ*/final_output/consensus.fasta | awk '/^>/ {split($0,a,"_"); print ">"a[2]} !/^>/ {print}' > genomes/$DT/tmp_raw.fasta

# Get headers from sequences.
grep ">" genomes/$DT/tmp_raw.fasta | sed 's/^>//' > genomes/$DT/tmp_header.txt

# Export only ones not exported before. 
#if [ -f export/exported.txt ]; then
#  comm -13 <(sort -u export/exported.txt) <(sort -u export/$DT/tmp_header.txt) > export/$DT/tmp_notexported.txt # Find LIB-IDs that have not been exported.
#else
  cat genomes/$DT/tmp_header.txt > genomes/$DT/tmp_notexported.txt
#  echo -e "library_id" > export/exported.txt
#fi
#  
#if [ ! -s export/$DT/tmp_notexported.txt ]; then
#  echo "All sequences has already been exported. Exiting."
#  rm -r export/$DT
#  exit 1
#fi

# Retain only lines from MAP_ID that are not exported.
awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if ($1 in seqs) {print $0}}' genomes/$DT/tmp_notexported.txt $MAP_ID > genomes/$DT/tmp_toexport.txt

# Add new name to "tmp_toexport.txt".
awk -F'\t' -v nmstring=$NAMESTRING '{newID=nmstring; sub(/@/,$2,newID); print $0"\t"newID }' genomes/$DT/tmp_toexport.txt > tmp && mv tmp genomes/$DT/tmp_toexport.txt

# Subset to sequences that should be exported and rename according to the new name.
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' genomes/$DT/tmp_raw.fasta | # tabularize.
awk -F'\t' 'FNR == NR {seqs[$1]=$2; next} {if (">"$1 in seqs) {print ">"$3"\t"seqs[">"$1]} else {print $0" had no matching sequence, excluding." > "/dev/stderr"}}' - genomes/$DT/tmp_toexport.txt > genomes/$DT/sequences.fasta # Subset to ones to export.

# Check for duplicates and remove with warning.
cut -f1 genomes/$DT/sequences.fasta | sort | uniq -c | awk '$1 > 1 {print $2}' > genomes/$DT/tmp_dups.txt

while read i; do
  echo "$i appears more than once, removing."
done < genomes/$DT/tmp_dups.txt

grep -Fv -f genomes/$DT/tmp_dups.txt genomes/$DT/sequences.fasta | # remove.
tr "\t" "\n" > tmp && mv tmp genomes/$DT/sequences.fasta # de-tabularize.

###############################################################################
# Dump metadata.
###############################################################################

echo -e "library_id\tssi_id\tnew_id" > genomes/$DT/idmapping.tsv
cat genomes/$DT/tmp_toexport.txt >> genomes/$DT/idmapping.tsv

# Clean-up folder.
rm genomes/$DT/tmp_*

exit 1

###############################################################################
# Remove human reads.
###############################################################################

if [ ! -z ${DUMPREADS+x} ]; then
  echo "Output .fastq file for each genome"
  
  # List all available files.
  for i in "${ARTICDIR[@]}"; do
    for j in $i/articminion/*.consensus.fasta; do echo $j; done
  done | sed 's/.consensus.fasta/.sorted.bam/' - > export/$DT/tmp_bamfiles
  
  # for each mapping, run sanitize-me.
  awk -v outdir=export/$DT -v human=$HUMANREF '{print $1":"outdir":"human}' export/$DT/tmp_toexport.txt | parallel -j $THREADS --colsep ':' --bar \
  '
  # Get the corresponding .bam file.
  BAMFILE=$(grep {1} {2}/tmp_bamfiles)
  
  if [ -z "$BAMFILE" ]; then
    >&2 echo "warning: .bam file for {1} was not in the artic output." 
  else
    # Extract mapped reads and remove human reads with the CDC protocol: https://github.com/CDCgov/SanitizeMe
    samtools fastq --threads 1 -F 4 $BAMFILE 2> /dev/null |\
    minimap2 -ax map-ont {3} - -t 1 2> /dev/null |\
    samtools view --threads 1 -u -f 4 - 2> /dev/null |\
    samtools bam2fq --threads 1 - 2> /dev/null |\
    gzip -c - > {2}/mapped_fastq/"$(sed s/\\//_/g <<< {1})".fastq.gz 2> /dev/null
  fi' 
  
  echo "md5 checksum"
  
  # md5 checksum.
  md5sum export/$DT/mapped_fastq/*fastq.gz | 
  awk '{n = split($2, a, "/"); print $1"\t"a[n]}' - | # remove file path.
  awk -F '\t' '{$3=$2; sub(/.fastq.gz/,"",$3); print $3"\t"$1"\t"$2}' - > export/$DT/tmp_md5sums.tsv # make ID column.  
  
  # Add to metadata.
  join <(sort genomes/$DT/tmp_toexport.txt) <(sort genomes/$DT/tmp_md5sums.tsv) -t $'\t' >> export/$DT/metadata.tsv
fi

# Write the parsed LIB-IDs to "exported.txt".
#cut -f 1 genomes/$DT/metadata.tsv | sed 1d - >> genomes/exported.txt


