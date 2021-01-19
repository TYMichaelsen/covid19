#!/usr/bin/env bash
# Written by Vang Le-Quy
# Date: 20200506
# A script to extract sequences from a subset of metadata downloaded from auspice

THISDIR=$(dirname $(readlink -f $0))
source ${THISDIR}/utils.sh
DISTDIR="/srv/rbd/covid19"
NEXTSTRAINOUT="${DISTDIR}/nextstrain"

SUBMETA=${1:-"none"}

GENOMES=${2:-$(findTheLatest "${DISTDIR}/genomes/*export")/sequences.fasta}
GENOUT=${3:-${SUBMETA/.tsv}.genomes.fasta}

if [ $SUBMETA = "none" ]; then
    echo "Subset of metadata is required"
    echo "Usage: $0 <subset_metadata> <genome_file> <output_subset_genome_file>"
    exit 1
fi

IMGDIR="/srv/rbd/thecontainer"
LASTESTIMG=$(findTheLatest "${IMGDIR}/*sif")
SINGIMG=${4:-$LASTESTIMG}

echo Using Singularity image: $SINGIMG
# DISTDIR="/srv/rbd/covid19/current"

echo "Script to extract genome sequences by using a subset of metadata"


# echo -e singularity exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "awk '(NR >1){print \$1}' $META|seqtk subseq $GENOMES - |seqtk seq -Cl60 >$GENOUT"
singularity exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "awk '(NR >1){print \$1}' $SUBMETA|seqtk subseq $GENOMES - |seqtk seq -Cl60 >$GENOUT"

echo Done extracting sequences to $GENOUT
