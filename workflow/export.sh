#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-s file -m file -a string -o dir -t int] 
-- COVID-19 data export for GISAID v. $VERSION:  

Arguments:
    -h  Show this help text.
    -m  Metadata .tsv file.
    -t  Number of threads.

Output:
    To come.

Note: 
    (DEPRECIATED) Using the option -H it is possible to specify an extended sample name that should be used instead of just the ID. 
    Example: 
    ID = HH-1
    -H = "DK_ALAB-@_2020" will lead to renaming as 'DK_ALAB-SSI-HH-1_2020'
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hs:m:a:o:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) META=$OPTARG;;
    t) THREADS=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${META+x} ]; then echo "-m $MISSING"; exit 1; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;

if [ ! -f $META ]; then echo "-m does not exist. It has likely been updated, check again."; exit 1; fi

### Code.----------------------------------------------------------------------
# Set dependent directories/data.
SEQS=QC/genomes/raw.fasta
ARTICDIR=( CJ* )
HUMANREF=/srv/rbd/covid19/data/reference/human_g1k_v37.fasta
NAMESTRING="DK/ALAB-@/2020"

# Make time-stamped export subfolder.
DT=$(date +%Y_%m_%d_%H-%M)"_export"
mkdir -p export/$DT
mkdir -p export/$DT/mapped_fastq

# Import metadata.
awk -F '\t' 'NR == 1 {for (i=1; i<=NF; i++) ix[$i] = i} NR > 1 {print $ix["ssi_id"]"\t"$ix["library_id"]"\t"$ix["curate_exclude"]"\t"$ix["sample_date"]}' $META | # Get the right columns.
awk -F '\t' '$2 != ""' - > export/$DT/tmp_all.txt
awk '$3 == "Yes" {print}' export/$DT/tmp_all.txt > export/$DT/tmp_exclude.txt
awk '$3 != "Yes" {print}' export/$DT/tmp_all.txt > export/$DT/tmp_pass.txt

# Export only ones not exported before. 
if [ -f export/exported.txt ]; then
  comm -13 <(sort -u <(cut -f2 export/exported.txt)) <(sort -u <(cut -f2 export/$DT/tmp_pass.txt)) > export/$DT/tmp_toexport.txt # Find LIB-IDs that have not been exported.
else
  cut -f2 export/$DT/tmp_pass.txt > export/$DT/tmp_toexport.txt
  echo -e "ssi_id\tlibrary_id\tfinal_id" > export/exported.txt
fi
  
if [ ! -s export/$DT/tmp_toexport.txt ]; then
  echo "All sequences specified in -m has already been exported. Exiting."
  exit 1
fi
 
EXPORTED=$(comm -13 <(sort -u export/$DT/tmp_toexport.txt) <(sort -u <(cut -f2 export/$DT/tmp_pass.txt)) | wc -l | sed 's/ .*//')
echo "$EXPORTED of the sequences in -m has been exported before." 
echo "Start exporting $(wc -l export/$DT/tmp_toexport.txt | sed 's/ .*//') of the $(wc -l export/$DT/tmp_pass.txt | sed 's/ .*//') sequences in -m passing manual QC."

# Subset to ones not exported before.
awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if ($2 in seqs) {print $0}}' export/$DT/tmp_toexport.txt export/$DT/tmp_pass.txt > tmp && mv tmp export/$DT/tmp_toexport.txt

# Subset seqs to metadata.
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' QC/genomes/raw.fasta | # tabularize.
awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if (">"$2 in seqs) {print seqs[">"$2]} else {print $0" had no matching sequence, excluding." > "/dev/stderr"}}' - export/$DT/tmp_toexport.txt | # Subset to ones in metadata.
tr "\t" "\n" > export/$DT/tmp_toexport.fasta # de-tabularize.

echo "Excluded $(wc -l export/$DT/tmp_exclude.txt | sed 's/ .*//') sequences."
echo "Continuing with $(grep ">" -c export/$DT/tmp_toexport.fasta) sequences."

# Add new name to "toexport.txt".
awk -F'\t' -v nmstring=$NAMESTRING '{newID=nmstring; sub(/@/,$1,newID); print $0"\t"newID }' export/$DT/tmp_toexport.txt > tmp && mv tmp export/$DT/tmp_toexport.txt

###############################################################################
# Remove human reads.
###############################################################################
echo "Output .fastq file for each genome"

# List all available files.
for i in "${ARTICDIR[@]}"; do
  for j in $i/articminion/*.consensus.fasta; do echo $j; done
done | sed 's/.consensus.fasta/.sorted.bam/' - > export/$DT/tmp_bamfiles

# for each mapping, run sanitizeme.
awk -v outdir=export/$DT -v human=$HUMANREF '{print $1":"$2":"$5":"outdir":"human}' export/$DT/tmp_toexport.txt | parallel -j $THREADS --colsep ':' --bar \
'
# Get the corresponding .bam file.
BAMFILE=$(grep {2} {4}/tmp_bamfiles)

if [ -z "$BAMFILE" ]; then
  >&2 echo "warning: .bam file for {1} (LIB-ID: {2}) was not in the artic output." 
else
  # Extract mapped reads and remove human reads with the CDC protocol: https://github.com/CDCgov/SanitizeMe
  samtools fastq --threads 1 -F 4 $BAMFILE 2> /dev/null |\
  minimap2 -ax map-ont {5} - -t 1 2> /dev/null |\
  samtools view --threads 1 -u -f 4 - 2> /dev/null |\
  samtools bam2fq --threads 1 - 2> /dev/null |\
  gzip -c - > {4}/mapped_fastq/"$(sed s/\\//_/g <<< {3})".fastq.gz 2> /dev/null
fi' 

echo "md5 checksum"

# md5 checksum.
md5sum export/$DT/mapped_fastq/*fastq.gz | 
awk '{n = split($2, a, "/"); print $1"\t"a[n]}' - | # remove file path.
awk -F '\t' '{$3=$2; sub(/.fastq.gz/,"",$3); gsub(/_/,"/",$3); print $3"\t"$1"\t"$2"\t"}' - > export/$DT/tmp_md5sums.tsv # make ID column.  

# Rename genomes.
awk -F'\t' 'FNR == NR{a[$2]=$5; next} {for (i in a) sub(i, a[i]); print}' export/$DT/tmp_toexport.txt export/$DT/tmp_toexport.fasta > export/$DT/sequences.fasta

###############################################################################
# Dump metadata.
###############################################################################

echo -e "ID\tsample_date\tmd5sum\tfastq_file" > export/$DT/metadata.tsv
join <(awk -F'\t' '{ print $5"\t"$4 }' export/$DT/tmp_toexport.txt | sort) <(sort export/$DT/tmp_md5sums.tsv) -t $'\t' >> export/$DT/metadata.tsv

# Write the passed ID, LIB-IDs, and new IDs to "exported.txt".
cut -f 1,2,5 export/$DT/tmp_toexport.txt >> export/exported.txt

# Clean-up folder.
mkdir -p export/$DT/tmpdir
mv export/$DT/tmp_* export/$DT/tmpdir

exit 1

###############################################################################
# Final cross-check and save which sequences were exported.
###############################################################################

# Get final IDs from all outputs.
ls -1 export/$DT/mapped_fastq | tr '\n' '\n' | sed 's/.fastq.gz//' | sed 's/_/\//g' > export/$DT/reads
grep ">" export/$DT/sequences.fasta | sed 's/^>//' > export/$DT/seqs
cut -f2 export/$DT/md5sums.txt | sed 's/.fastq.gz//' | sed 's/_/\//g' > export/$DT/md5sums

# Match them up and report missing. 
awk '{a[$1]++} END {for (i in a) print i}' export/$DT/reads export/$DT/seqs export/$DT/md5sums > export/$DT/union

awk 'FNR == NR {a[$1]=$1; next} { if ($1 in a) {print $1} else {print ""}' export/$DT/union export/$DT/reads


awk 'FNR==1{x++} {a[$1][x]} END {for(i in a){q=0;for(j in a[i]){q++};print q,q"/"x,b[i]}}' file* 

comm -13 <(sort -u <(grep ">" export/$DT/sequences.fasta | sed 's/^>//')) <(ls export/$DT/mapped_fastq/* | sed 's/_/\\//g')) > export/$DT/toexport.txt # Find LIB-IDs that have not been exported.

comm 




