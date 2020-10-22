#!/bin/bash 

# This script checks for new batches in /srv/rbd/covid19/processing. If a new batch is present, it monitors for a while - if the folder is static and has >20Gb it starts the covid19_workflow.sh script.

# List all batches processed so far. Criteria for proccesed is consensus.fasta in /final_output/ folder.
for batch in $(ls -dtr /srv/rbd/covid19/processing/?J*); do 
  if [ -f $batch/final_output/consensus.fasta ]; then
    echo $batch
  fi
done > /srv/rbd/covid19/processing/processed.txt

while : ; do

# Check if a new batch has arrived.
comm -23 <(ls -d /srv/rbd/covid19/processing/?J* | sort) <(sort /srv/rbd/covid19/processing/processed.txt) > /srv/rbd/covid19/processing/missing.txt

while [ -s /srv/rbd/covid19/processing/missing.txt ]; do 
  
  DIR=$(head -n 1 /srv/rbd/covid19/processing/missing.txt)
  
  echo "[$(date +"%b %d %T")] Found a new batch: $(basename $DIR)"
  
  # Check if the folder is static (i.e. no more output recived from elsewhere).
  
  while : ; do
    old_files=$(find $DIR)
    
    echo "[$(date +"%b %d %T")] Wait to check if $(basename $DIR) is static. Check again in 5 minutes."
    
    sleep 5m
    
    new_files=$(find $DIR)
    
    diff_files=$(diff -q <(echo $old_files) <(echo $new_files))
    
    if [ -z "$diff_files" ]; then 
      echo "[$(date +"%b %d %T")] $(basename $DIR) is static, continue with workflow processing." 
      break
    fi
  done  

  # Test if there is >20 Gb data. If not, cancel the job. 
  nfastq=$(find $DIR/fastq -name *fastq | wc -l)
  nbases=$(( $nfastq * 4000 * 1050 ))

  if [ $nbases -ge 20000000000 ]; then

    # Run the workflow.  
    covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1

    # Remove from top of missing.txt and add to processed.txt
    tail -n +2 /srv/rbd/covid19/processing/missing.txt > tmp && mv tmp /srv/rbd/covid19/processing/missing.txt
    cat /srv/rbd/covid19/processing/processed.txt <(echo $DIR) > tmp && mv tmp /srv/rbd/covid19/processing/processed.txt

  else
    
    echo "[$(date +"%b %d %T")] $(basename $DIR) is static but has <20 Gb data, adding it to bottom of missing.txt and continue." 
    
    # Remove from top of missing.txt and add to bottom of missing.txt
    tail -n +2 /srv/rbd/covid19/processing/missing.txt       > tmp && mv tmp /srv/rbd/covid19/processing/missing.txt
    cat /srv/rbd/covid19/processing/missing.txt <(echo $DIR) > tmp && mv tmp /srv/rbd/covid19/processing/missing.txt
  fi

done

echo "[$(date +"%b %d %T")] No new batches found. Check again in 30 minutes." 

sleep 30m

done