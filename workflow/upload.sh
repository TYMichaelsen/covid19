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
  Depends on finding the columns 'gisaid_id', 'library_id','curate_exclude','region', and 'sample_date' in the metadata. 
  These should be formatted as:
    gisaid_id      = id associated with the sample.
    library_id     = id associated with the library.
    curate_exclude = column indicating if the sample should be excluded. Empty entries means inclusion. 
    region         = column with sample DK region.
    sample_date    = column with sample date in YYYY-MM-DD format.

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

#META=/srv/rbd/covid19/upload/2020-07-16-08-49_upload/tmp_metadata.tsv
#SEQS=/srv/rbd/covid19/genomes/2020-07-16-08-49_export/sequences.fasta
#WORKFLOW_PATH=/srv/rbd/covid19/git/covid19/workflow

WORKFLOW_PATH="$(dirname "$(readlink -f "$0")")"
OUTDIR=/srv/rbd/covid19/upload/$DT

# Set dependent directories/data.
ARTICDIR=( /srv/rbd/covid19/processing/CJ*/processing )
HUMANREF=$WORKFLOW_PATH/dependencies/human_g1k_v37.fasta

# Make time-stamped export subfolder.
mkdir -p $OUTDIR
mkdir -p $OUTDIR/mapped_fastq

# Reorder columns in metadata.
awk -F '\t' 'NR == 1 {for (i=1; i<=NF; i++) ix[$i] = i} NR > 1 {print $ix["ssi_id"]"\t"$ix["library_id"]"\t"$ix["gisaid_id"]"\t"$ix["region"]"\t"$ix["sample_date"]}' $META > tmp && mv tmp $META

# Subset seqs to metadata.
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' $SEQS | # tabularize.
awk -F'\t' 'FNR == NR {seqs[$1]=$2; next} {if (">"$1 in seqs) {print ">"$3"\t"seqs[">"$1]} else {print $1" had no matching sequence, excluding." > "/dev/stderr"}}' - $META | # Subset to ones in metadata.
tr "\t" "\n" > $OUTDIR/sequences.fasta # de-tabularize.

echo "Continuing with $(grep ">" -c $OUTDIR/sequences.fasta) sequences."

###############################################################################
# Remove human reads.
###############################################################################
echo "Output .fastq file for each genome"

# List all available files.
for i in "${ARTICDIR[@]}"; do
  for j in $i/articminion/*.consensus.fasta; do echo $j; done
done | sed 's/.consensus.fasta/.sorted.bam/' - > $OUTDIR/tmp_bamfiles

head -n 10 $OUTDIR/tmp_bamfiles > tmp && mv tmp $OUTDIR/tmp_bamfiles

# Make input file to parallel.
awk -F'\t' -v outdir=$OUTDIR -v human=$HUMANREF '{print $3":"$2":"outdir":"human}' $OUTDIR/tmp_metadata.tsv > $OUTDIR/tmp_toparallel.txt

# Run in singularity.
SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
singularity --silent exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "OUTDIR=$OUTDIR; THREADS=$THREADS; . $WORKFLOW_PATH/parallel.snipet.sh"

echo "md5 checksum"

# md5 checksum.
md5sum $OUTDIR/mapped_fastq/*fastq.gz | 
awk '{n = split($2, a, "/"); print $1"\t"a[n]}' - | # remove file path.
awk -F '\t' '{$3=$2; sub(/.fastq.gz/,"",$3); gsub(/_/,"/",$3); print $3"\t"$1"\t"$2"\t"}' - > $OUTDIR/tmp_md5sums.tsv # make ID column.  

###############################################################################
# Dump metadata.
###############################################################################

echo -e "gisaid_id\tsample_date\tmd5sum\tfastq_file" > $OUTDIR/metadata.tsv
join <(awk -F'\t' '{ print $3"\t"$4 }' $OUTDIR/tmp_metadata.tsv | sort) <(sort $OUTDIR/tmp_md5sums.tsv) -t $'\t' >> $OUTDIR/metadata.tsv

# Clean-up folder.
rm $OUTDIR/tmp_* 
