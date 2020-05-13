#!/usr/bin/env bash
# A wrapper script to start Auspice server for Nextstrain
# Written by Vang Le-Quy
# Date 20200506


if [ x"$1" ==  "x-h" -o x"$1" == "x--help" ]; then
    echo "Usage: $0 [default|/path/to/auspice/output] [PORT]"
    echo "PORT is optional and random by default between 4001-6500"
    echo "If you specify PORT, it will be used for auspice"
    exit 1
fi


SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
MYPORT=$(shuf -i 4001-6500 -n 1)
MYPORT=${2:-$MYPORT}

THISDIR=$(dirname $(readlink -f $0))

# Source utility functions
source ${THISDIR}/utils.sh

LATESTOUTPUT=$(findTheLatest "/srv/rbd/covid19/nextstrain/*nextstrain")

DATDIR="${1:-${LATESTOUTPUT}/auspice}"

if [ "x${1}" == "xdefault" ]; then
    DATDIR="${LATESTOUTPUT}"
fi

DATDIR=$(readlink -f $DATDIR)

singularity exec -B $DATDIR:$DATDIR -B /srv/rbd:/srv/rbd $SINGIMG bash -c "source activate nextstrain; PORT=$MYPORT auspice view --datasetDir $DATDIR" 
