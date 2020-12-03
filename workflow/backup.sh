#!/bin/bash 

# This script takes a batch dir as input and then 1) make a backup to /srv/backup without rawdata and 2) remove rawdata from the batch dir.
# Currently only TYM has r/w permission to /srv/backup so he is the only one who can run the script.

# List all batches processed so far. Criteria for backed up proccesed is a non-empty consensus.fasta in /final_output/ folder.
batches=$(for batch in $(ls -dtr /srv/rbd/covid19/processing/?J*) $(ls -dtr /srv/data_1/?J*); do
  if [ -s $batch/final_output/consensus.fasta ]; then
    echo $batch
  fi
done) 

while : ; do

  # Go through each batch and rsync.
  for i in $batches; do 
    rsync -arv --log-file=/srv/backup/backup.log --exclude /rawdata --exclude /demultiplexed --exclude /filtered $i/ /srv/backup/$(basename $i)
  done 
  
  # Wait 2 hours and do it again.
  sleep 2h

done
  