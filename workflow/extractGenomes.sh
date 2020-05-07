#!/usr/bin/env bash
# Written by Vang Le-Quy
# Date: 20200506
# A script to extract sequences from a subset of metadata downloaded from auspice

SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
SUBMETA="$1"
NEXTSTRAIN="/srv/rbd/covid19/nextstrain/2020-05-05_nextstrain_clean"
GENOMES=${2:-$NEXTSTRAIN/raw.fasta}

echo "Script to extract genome sequences by using a subset of metadata downloaded from Auspice"

if [ $# -lt 1 ]; then
    echo "Usage: $0 <submeta.tsv> [optional: genomes.fasta]"
    echo "By default, $NEXTSTRAIN/raw.fasta is used to find genomes"
    exit 1
fi

GENOUT=${SUBMETA/.tsv}.genomes.fasta

singularity exec -B /srv/rbd:/srv/rbd $SINGIMG bash -c "awk '(NR >1){print \$1}' $SUBMETA|seqtk subseq $GENOMES - |seqtk seq -Cl60 >$GENOUT"

echo Done extracting sequences to $GENOUT
