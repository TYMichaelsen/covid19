#!/bin/bash
# By Thomas Y. Michaelsen
VERSION=0.1.0

### Description ----------------------------------------------------------------

USAGE="$(basename "$0") [-h] [-s file -m file -t int] 
-- COVID-19 data export for GISAID v. $VERSION:  

Arguments:
    -h  Show this help text.
    -m  Metadata .tsv file.
    -t  Optional: if provided the name of folder. 

Output:
    To come.

Note: 
    The script extracts the mapping between library_id and user-specified id. Samples not in library_id will be excluded. It also checks for duplicates and remove these.
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hm:t:' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    m) META=$OPTARG;;
    t) DT=$OPTARG;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${META+x} ]; then echo "-m $MISSING"; exit 1; fi;
if [ -z ${DT+x} ]; then DT=$(date +%Y_%m_%d_%H-%M)"_export"; fi;

### Code.----------------------------------------------------------------------

GITDIR="$(dirname "$(readlink -f "$0")")"

# How the header should look like.
NAMESTRING="hCoV-19/DK/ALAB-@/2020"

# Make time-stamped export subfolder.
mkdir -p genomes
mkdir -p genomes/$DT

# Load and parse metadata.
awk -F '\t' 'NR == 1 {for (i=1; i<=NF; i++) ix[$i] = i} NR > 1 {print $ix["library_id"]"\t"$ix["ssi_id"]}' $META | # Get the right columns.
awk -F '\t' '$1 != ""' - | # Keep only ones with library_id.
sed '/POS/d' - | # Remove all controls.
sed '/NEG/d' - > genomes/$DT/tmp_mapping.tsv # Remove all controls.

# Sweep batch-folders for consensus sequences (OBS: renaming headers according to AAU library naming scheme).
cat processing/?J*/final_output/consensus.fasta | awk '/^>/ {split($0,a,"_"); print ">"a[2]} !/^>/ {print}' > genomes/$DT/tmp_raw.fasta

# Get headers from sequences.
grep ">" genomes/$DT/tmp_raw.fasta | sed 's/^>//' > genomes/$DT/tmp_header.txt

cat genomes/$DT/tmp_header.txt > genomes/$DT/tmp_notexported.txt

# Export only lines from tmp_mapping.tsv.
awk -F'\t' 'FNR == NR {seqs[$1]=$0; next} {if ($1 in seqs) {print $0}}' genomes/$DT/tmp_notexported.txt genomes/$DT/tmp_mapping.tsv > genomes/$DT/tmp_toexport.txt

# Add new name to "tmp_toexport.txt".
awk -F'\t' -v nmstring=$NAMESTRING '{newID=nmstring; sub(/@/,$2,newID); print $0"\t"newID }' genomes/$DT/tmp_toexport.txt > tmp && mv tmp genomes/$DT/tmp_toexport.txt

# Subset to sequences that should be exported and rename according to the new name.
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' genomes/$DT/tmp_raw.fasta | # tabularize.
awk -F'\t' 'FNR == NR {seqs[$1]=$2; next} {if (">"$1 in seqs) {print ">"$3"\t"seqs[">"$1]} else {print $0" had no matching sequence, excluding." > "/dev/stderr"}}' - genomes/$DT/tmp_toexport.txt > genomes/$DT/sequences.fasta # Subset to ones to export.

# Check for duplicates and remove with warning.
cut -f1 genomes/$DT/sequences.fasta | sort | uniq -c | awk '$1 > 1 {sub(/^>/,"",$2); print $2}' > genomes/$DT/tmp_dups.txt

while read i; do
  echo "$i appears more than once, removing."
done < genomes/$DT/tmp_dups.txt

grep -Fv -f genomes/$DT/tmp_dups.txt genomes/$DT/sequences.fasta | # remove.
tr "\t" "\n" > tmp && mv tmp genomes/$DT/sequences.fasta # de-tabularize.

###############################################################################
# Dump metadata.
###############################################################################

echo -e "library_id\tssi_id\tgisaid_id" > genomes/$DT/idmapping.tsv
cat genomes/$DT/tmp_toexport.txt >> genomes/$DT/idmapping.tsv

# Also remove duplicates from id-mapping.
grep -Fv -f genomes/$DT/tmp_dups.txt genomes/$DT/idmapping.tsv > tmp && mv tmp genomes/$DT/idmapping.tsv

# Clean-up folder.
rm genomes/$DT/tmp_*
