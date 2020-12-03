#!/bin/bash 

# Need to work here to get permission to make temporary files.
cd /srv/rbd/covid19/processing

# This script checks for new batches in:
# - /srv/rbd/covid19/processing
# - /srv/data_1
# - /srv/data_1/gridion

# If a new batch is present, it monitors for a while -
# if the folder is static and has enough data it starts the covid19_workflow.sh script.

# Logging
LOG_NAME="/srv/rbd/covid19/processing/autorun_log_$(date +"%Y-%m-%d_%H-%M").txt"
echo "autorun log" >> $LOG_NAME
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

# List all batches processed so far. Criteria for proccesed is a non-empty consensus.fasta in /final_output/ folder.
for batch in $(ls -dtr /srv/rbd/covid19/processing/?J*) $(ls -dtr /srv/data_1/?J*) $(ls -dtr /srv/data_1/gridion/?J*); do
  if [ -s $batch/final_output/consensus.fasta ]; then
    echo $batch
  fi
done > /srv/rbd/covid19/processing/processed.txt

chmod 777 /srv/rbd/covid19/processing/processed.txt

while : ; do

# Check if a new batch has arrived.
comm -23 <(ls -d /srv/rbd/covid19/processing/?J* /srv/data_1/?J* /srv/data_1/gridion/?J* | sort) <(sort /srv/rbd/covid19/processing/processed.txt) > /srv/rbd/covid19/processing/missing.txt

chmod 777 /srv/rbd/covid19/processing/missing.txt

while [ -s /srv/rbd/covid19/processing/missing.txt ]; do 
  
  DIR=$(head -n 1 /srv/rbd/covid19/processing/missing.txt)
  
  echo "[$(date +"%b %d %T")] Found a new batch: $(basename $DIR)"
  
  ### Check if the folder is static (i.e. no more output recived from elsewhere).------------------------------------
  old_files=$(find $DIR)
    
  echo "[$(date +"%b %d %T")] Wait 2 min to check if $(basename $DIR) is static."
    
  sleep 2m
    
  new_files=$(find $DIR)
   
  diff_files=$(diff -q <(echo $old_files) <(echo $new_files))
    
  if [ -z "$diff_files" ]; then 
    echo "[$(date +"%b %d %T")] $(basename $DIR) is static, continue with workflow processing." 

    ### Test if there is enough data. If not, cancel the job.----------------------------------------------------------- 
    nfastq=$(find $DIR/fastq -name *fastq | wc -l)
    nbases=$(( $nfastq * 4000 * 1050 ))
    
    # Grep the number Gb specified by user. If none, the default 20Gb are used.
    if [ -s $DIR/datacap.txt ]; then DATACAP=$(($(cat $DIR/datacap.txt) * 1000000000)); else DATACAP=20000000000; fi
    
    if [ $nbases -ge $DATACAP ]; then

      # Run the workflow. 
      cd $(dirname $DIR) 
      covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1

      # Remove from top of missing.txt and add to processed.txt
      tail -n +2 /srv/rbd/covid19/processing/missing.txt > tmp_auto && mv tmp_auto /srv/rbd/covid19/processing/missing.txt
      cat /srv/rbd/covid19/processing/processed.txt <(echo $DIR) > tmp_auto && mv tmp_auto /srv/rbd/covid19/processing/processed.txt

    else
    
      echo "[$(date +"%b %d %T")] $(basename $DIR) is static but has <20 Gb data, adding it to bottom of missing.txt and continue." 
    
      # Remove from top of missing.txt and add to bottom of missing.txt
      tail -n +2 /srv/rbd/covid19/processing/missing.txt       > tmp_auto && mv tmp_auto /srv/rbd/covid19/processing/missing.txt
      cat /srv/rbd/covid19/processing/missing.txt <(echo $DIR) > tmp_auto && mv tmp_auto /srv/rbd/covid19/processing/missing.txt
    fi
  else 
    echo "[$(date +"%b %d %T")] $(basename $DIR) was not static, adding it to bottom of missing.txt and continue." 

    # Remove from top of missing.txt and add to bottom of missing.txt
    tail -n +2 /srv/rbd/covid19/processing/missing.txt       > tmp_auto && mv tmp_auto /srv/rbd/covid19/processing/missing.txt
    cat /srv/rbd/covid19/processing/missing.txt <(echo $DIR) > tmp_auto && mv tmp_auto /srv/rbd/covid19/processing/missing.txt
  fi

done

echo "[$(date +"%b %d %T")] No new batches found. Check again in 30 minutes." 

sleep 30m

done