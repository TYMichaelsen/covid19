#!/usr/bin/env bash

echo "Example run:"
echo "$0 /srv/rbd/covid19/current/export/2020_04_28_12-27_export/sequences.fasta /srv/rbd/covid19/current/metadata/metadata_SSI/2020-04-28_metadata_nextstrain.tsv  ./TestNextStrain"


SINGIMG="/srv/rbd/covid19/thecontainer/covid19_latest.sif"
# DISTDIR="/srv/rbd/covid19/current"
DISTDIR="/srv/rbd"
GENOMEFASTA="$1"
METAFILE="$2"
OUTDIR="${3:-$PWD/nextstrainRun}"
THISDIR=$(dirname $(readlink -f $0))
NEXTSTRAIN_SCRIPT=${THISDIR}/nextstrain.sh

echo "Running comand:"
echo "singularity exec  -B $DISTDIR:$DISTDIR -B $HOME:$HOME $SINGIMG bash -c \"source activate nextstrain; $NEXTSTRAIN_SCRIPT -s $GENOMEFASTA -m $METAFILE -o $OUTDIR\""
singularity exec  -B $DISTDIR:$DISTDIR -B $HOME:$HOME $SINGIMG bash -c "source activate nextstrain; $NEXTSTRAIN_SCRIPT -s $GENOMEFASTA -m $METAFILE -o $OUTDIR"
