#!/usr/bin/env bash
SINGIMG="/srv/rbd/covid19/thecontainer/covid19_latest.sif"
MYPORT=$(shuf -i 4000-6500 -n 1)
DATDIR=${1:-/srv/rbd/covid19/current/analysis/nextstrain/auspice/}
singularity exec -B $DATDIR:$DATDIR -B /srv/rbd:/srv/rbd $SINGIMG bash -c "source activate nextstrain; PORT=$MYPORT auspice view --datasetDir $DATDIR"
