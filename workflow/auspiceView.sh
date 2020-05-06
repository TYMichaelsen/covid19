#!/usr/bin/env bash
# A wrapper script to start Auspice server for Nextstrain
# Written by Vang Le-Quy
# Date 20200506


if [ $# -lt 1 ]; then
    echo "Usage: $0 [default|/path/to/auspice/output] [PORT]"
    echo "PORT is optional and random by default between 4001-6500"
    echo "If you specify PORT, it will be used for auspice"
    exit 1
fi


SINGIMG="/srv/rbd/thecontainer/covid19_latest.sif"
MYPORT=$(shuf -i 4001-6500 -n 1)
MYPORT=${2:-$MYPORT}

DATDIR="${1:-/srv/rbd/covid19/nextstrain/2020-05-05_nextstrain_clean/auspice}"

if [ "x${1}" == "xdefault" ]; then
    DATDIR="/srv/rbd/covid19/nextstrain/2020-05-05_nextstrain_clean/auspice"
fi

DATDIR=$(readlink -f $DATDIR)

singularity exec -B $DATDIR:$DATDIR -B /srv/rbd:/srv/rbd $SINGIMG bash -c "source activate nextstrain; PORT=$MYPORT auspice view --datasetDir $DATDIR" 
