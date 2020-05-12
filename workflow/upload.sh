#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-s file -m file -a string -o dir -t int] 
-- COVID-19 prep data for upload v. $VERSION:  

Arguments:
    -h  Show this help text.
    -m  Metadata .tsv file.
    -s  Sequences.
    -t  Number of threads.
    -x  Optional: if set, the name of output dir.

Description:
  Dumps timestamped folder in upload/ containing sequences, fastq reads, and metadata.
  Depends on finding the columns 'gisaid_id', 'library_id','curate_exclude', and 'SampleDate' in the metadata. 
  These should be formatted as:
    gisaid_id      = id associated with the sample.
    library_id     = id associated with the library.
    curate_exclude = column indicating if the sample should be excluded. Empty entries means inclusion. 
    SampleDate     = column with sample date in YYYY-MM-DD format.

"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hm:s:t:x:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) META=$OPTARG;;
    s) SEQS=$OPTARG;;
    t) THREADS=$OPTARG;;
    x) DT=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${META+x} ]; then echo "-m $MISSING"; exit 1; fi;
if [ -z ${SEQS+x} ]; then echo "-s $MISSING"; exit 1; fi;
if [ -z ${THREADS+x} ]; then THREADS=50; fi;
if [ -z ${DT+x} ]; then DT=$(date +%Y-%m-%d-%H-%M)"_upload"; fi;

if [ ! -f $META ]; then echo "-m does not exist. It has likely been updated, check again."; exit 1; fi
if [ ! -f $SEQS ]; then echo "-s does not exist."; exit 1; fi

### Code.----------------------------------------------------------------------

WORKFLOW_PATH="$(dirname "$(readlink -f "$0")")"

# Set dependent directories/data.
ARTICDIR=( processing/CJ*/processing )
HUMANREF=$WORKFLOW_PATH/dependencies/human_g1k_v37.fasta

# Make time-stamped export subfolder.
mkdir -p upload/$DT
mkdir -p upload/$DT/mapped_fastq

# Import metadata.
awk -F '\t' 'NR == 1 {for (i=1; i<=NF; i++) ix[$i] = i} NR > 1 {print $ix["gisaid_id"]"\t"$ix["library_id"]"\t"$ix["curate_exclude"]"\t"$ix["SampleDate"]}' $META | # Get the right columns.
awk -F '\t' '$2 != "" && $1 != ""' - > upload/$DT/tmp_all.txt
awk -F '\t' '$3 != "" {print}' upload/$DT/tmp_all.txt > upload/$DT/tmp_exclude.txt
awk -F '\t' '$3 == "" {print}' upload/$DT/tmp_all.txt > upload/$DT/tmp_pass.txt

# Export only ones not exported before. 
if [ -f upload/exported.txt ]; then
  comm -13 <(sort -u <(cut -f2 upload/exported.txt)) <(sort -u <(cut -f2 upload/$DT/tmp_pass.txt)) > upload/$DT/tmp_toexport.txt # Find LIB-IDs that have not been exported.
else
  cut -f2 upload/$DT/tmp_pass.txt > upload/$DT/tmp_toexport.txt
  echo -e "ssi_id\tlibrary_id\tfinal_id" > upload/exported.txt
fi
  
if [ ! -s upload/$DT/tmp_toexport.txt ]; then
  echo "All sequences specified in -m has already been exported. Exiting."
  exit 1
fi

# Print out what to export.
EXPORTED=$(comm -13 <(sort -u upload/$DT/tmp_toexport.txt) <(sort -u <(cut -f2 upload/$DT/tmp_pass.txt)) | wc -l | sed 's/ .*//')
echo "$EXPORTED of the sequences in -m has been exported before." 
echo "Start exporting $(wc -l upload/$DT/tmp_toexport.txt | sed 's/ .*//') of the $(wc -l upload/$DT/tmp_pass.txt | sed 's/ .*//') sequences in -m passing manual QC."

# Subset to ones not exported before.
awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if ($2 in seqs) {print $0}}' upload/$DT/tmp_toexport.txt upload/$DT/tmp_pass.txt > tmp && mv tmp upload/$DT/tmp_toexport.txt

# Subset seqs to metadata.
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' $SEQS | # tabularize.
awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if (">"$1 in seqs) {print seqs[">"$1]} else {print $0" had no matching sequence, excluding." > "/dev/stderr"}}' - upload/$DT/tmp_toexport.txt | # Subset to ones in metadata.
tr "\t" "\n" > upload/$DT/sequences.fasta # de-tabularize.

echo "Continuing with $(grep ">" -c upload/$DT/sequences.fasta) sequences."

###############################################################################
# Remove human reads.
###############################################################################
echo "Output .fastq file for each genome"

# List all available files.
for i in "${ARTICDIR[@]}"; do
  for j in $i/articminion/*.consensus.fasta; do echo $j; done
done | sed 's/.consensus.fasta/.sorted.bam/' - > upload/$DT/tmp_bamfiles

# Make input file to parallel.
awk -F'\t' -v outdir=upload/$DT -v human=$HUMANREF '{print $1":"$2":"outdir":"human}' upload/$DT/tmp_toexport.txt > upload/$DT/tmp_toparallel.txt

# Run in singularity.
SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
singularity --silent exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "DT=$DT; THREADS=$THREADS; . $WORKFLOW_PATH/parallel.snipet.sh"

echo "md5 checksum"

# md5 checksum.
md5sum upload/$DT/mapped_fastq/*fastq.gz | 
awk '{n = split($2, a, "/"); print $1"\t"a[n]}' - | # remove file path.
awk -F '\t' '{$3=$2; sub(/.fastq.gz/,"",$3); gsub(/_/,"/",$3); print $3"\t"$1"\t"$2"\t"}' - > upload/$DT/tmp_md5sums.tsv # make ID column.  

###############################################################################
# Dump metadata.
###############################################################################

echo -e "gisaid_id\tsample_date\tmd5sum\tfastq_file" > upload/$DT/metadata.tsv
join <(awk -F'\t' '{ print $1"\t"$4 }' upload/$DT/tmp_toexport.txt | sort) <(sort upload/$DT/tmp_md5sums.tsv) -t $'\t' >> upload/$DT/metadata.tsv

# Write the passed gisaid_id and library_id to "exported.txt".
cut -f 1,2 upload/$DT/tmp_toexport.txt >> upload/exported.txt

# Clean-up folder.
rm upload/$DT/tmp_* 
