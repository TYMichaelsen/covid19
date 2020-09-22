#!/bin/bash 

# This script checks for new batches in /srv/rbd/covid19/processing. If a new batch is present, it monitors for a while - if the folder is static it starts the covid19_workflow.sh script.
#
# 

# FIRST USE ONLY: Make a date-sorted text file with all the batches processed so far.
# ls -dtr /srv/rbd/covid19/processing/?J* > /srv/rbd/covid19/processing/processed.txt

while : ; do

# Check if a new batch has arrived.
comm -23 <(ls -d /srv/rbd/covid19/processing/?J* | sort) <(sort /srv/rbd/covid19/processing/processed.txt) > /srv/rbd/covid19/processing/missing.txt

while [ -s /srv/rbd/covid19/processing/missing.txt ]; do 
  
  DIR=$(head -n 1 /srv/rbd/covid19/processing/missing.txt)
  
  echo "[$(date +"%T")] Found a new batch: $(basename $DIR)"
  
  # Check if the folder is static (i.e. no more output recived from AI node).
  
  while : ; do
    old_files=$(find $DIR)
    
    echo "[$(date +"%T")] Wait to check if $(basename $DIR) is static. Check again in 5 minutes."
    
    sleep 2
    
    new_files=$(find $DIR)
    
    diff_files=$(diff -q <(echo $old_files) <(echo $new_files))
    
    if [ -z "$diff_files" ]; then 
      echo "[$(date +"%T")] $(basename $DIR) is static, continue with workflow processing." 
      break
    fi
  done
  
  # Run the workflow.  
  covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1
  
  # Remove from missing.txt and add to processed.txt
  tail -n +2 /srv/rbd/covid19/processing/missing.txt > /srv/rbd/covid19/processing/missing.txt
  cat /srv/rbd/covid19/processing/processed.txt <(echo $DIR) > tmp && mv tmp /srv/rbd/covid19/processing/processed.txt
   
done

echo "[$(date +"%T")] No new batches found. Check again in 30 minutes." 

sleep 30m

done