#!/bin/bash 

# Need to work here to get permission to make temporary files.
cd /srv/rbd/covid19/processing

# This script checks for new batches in:
# - /srv/rbd/covid19/processing
# - /srv/data_1
# - /srv/data_1/gridion

# If a new batch is present, it monitors for a while -
# if the folder has enough data it starts the covid19_workflow.sh script.

# How many minutes between echo'ing that the process is still going.
LIFESIGN=30 
ATTEMPT=0

rm missing.txt
rm processed.txt

# Logging
LOG_NAME="/srv/rbd/covid19/processing/autorun_log_$(date +"%Y-%m-%d_%H-%M").txt"
echo "autorun log" >> $LOG_NAME
echo "Command: $0 $*" >> $LOG_NAME
exec &> >(tee -a "$LOG_NAME")
exec 2>&1

chmod 777 $LOG_NAME

echo "[$(date +"%b %d %T")] Started waiting for new batches. Will check every 2 min." 

while : ; do

# List all batches processed so far. Criteria for proccesed is a non-empty consensus.fasta in /final_output/ folder.
for batch in $(ls -dtr /srv/rbd/covid19/processing/?J*/) $(ls -dtr /srv/data_1/?J*/) $(ls -dtr /srv/data_1/gridion/?J*/); do
  if [ -s $batch/final_output/consensus.fasta ]; then
    echo $batch
  fi
done > /srv/rbd/covid19/processing/processed.txt

chmod 777 /srv/rbd/covid19/processing/processed.txt

# List all batches NOT processed.
comm -23 <(ls -d /srv/rbd/covid19/processing/?J*/ /srv/data_1/?J*/ /srv/data_1/gridion/?J*/ | sort) <(sort /srv/rbd/covid19/processing/processed.txt) > /srv/rbd/covid19/processing/missing.txt

chmod 777 /srv/rbd/covid19/processing/missing.txt

if [ -s /srv/rbd/covid19/processing/missing.txt ]; then 

  ATTEMPT=0

  while read DIR; do 
  
    echo "[$(date +"%b %d %T")] Found a new batch: $(basename $DIR)"
    
    ### Check if the folder is static (i.e. no more output recived from elsewhere).------------------------------------
    old_files=$(find $DIR)
    
    echo "[$(date +"%b %d %T")] Wait 2 min to check if $(basename $DIR) is static."

    sleep 2m
    
    new_files=$(find $DIR)
   
    diff_files=$(diff -q <(echo $old_files) <(echo $new_files))
    
    if [ -z "$diff_files" ]; then 

      ### Test if there is enough data. If not, cancel the job.----------------------------------------------------------- 
      nfastq=$(find $DIR/fastq -name *fastq | wc -l)
      nbases=$(( $nfastq * 4000 * 1050 ))
    
      # Grep the number Gb specified by user. If none, the default 20Gb are used.
      if [ -s $DIR/datacap.txt ]; then DATACAP=$(($(cat $DIR/datacap.txt) * 1000000000)); else DATACAP=20000000000; fi
    
      if [ $nbases -ge $DATACAP ]; then
        
        if [ -z "$(pgrep -fl covid19_workflow.sh)" ]; then
          
          # Run the workflow. 
          echo "[$(date +"%b %d %T")] $(basename $DIR) is static, continue with workflow processing." 
          cd $(dirname $DIR) 
          covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1 >/dev/null 2>&1
          echo "[$(date +"%b %d %T")] workflow completed for $(basename $DIR)." 
        else 
          echo "[$(date +"%b %d %T")] covid19_workflow.sh is running elsewhere. Skipping and continue." 
        fi  
      else
        echo "[$(date +"%b %d %T")] $(basename $DIR) is static but has <20 Gb data, skipping and continue." 
      fi
    else 
      echo "[$(date +"%b %d %T")] $(basename $DIR) was not static, skipping and continue." 
    fi
  
  done < /srv/rbd/covid19/processing/missing.txt

  echo "[$(date +"%b %d %T")] Started waiting for new batches. Will check every 2 min." 

else 
  ATTEMPT=$(($ATTEMPT+1))
  
  if [ ! $ATTEMPT -lt $LIFESIGN ]; then
    echo "[$(date +"%b %d %T")] Still waiting!"
    ATTEMPT=0
  fi
  sleep 1m
fi

done