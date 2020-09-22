#!/bin/bash 

# This script takes a batch dir as input and then 1) make a backup to /srv/backup without rawdata and 2) remove rawdata from the batch dir.
# Currently only TYM has r/w permission to /srv/backup so he is the only one who can run the script.

batch=$1

# Remove /rawdata/

rm -r $batch/rawdata

# rsync to /srv/backup

rsync -arv $batch /srv/backup/$(basename $batch)

# Print out batches that has not been backuped yet.

echo "Following batches has not been backed up yet:"
comm -23 <(ls -d /srv/rbd/covid19/processing/?J* | sort) <(ls -d /srv/backup/?J* | sort)
