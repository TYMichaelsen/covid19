#!/bin/bash 
FORCE_START=${1:-2500}
LAST_TIME=$(($(date +%H%M) - 2 ))

# Need to work here to get permission to make temporary files.
cd /srv/rbd/covid19/processing

# This script checks for new batches in:
# - /srv/rbd/covid19/processing
# - /srv/data_1
# - /srv/data_1/gridion

# If a new batch is present, it monitors for a while -
# if the folder has enough data it starts the covid19_workflow.sh script.

# How many minutes between echo'ing that the process is still going.
LIFESIGN=1
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
  
  CUR_TIME=$(date +%H%M)

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

    # Before doing anything, check that timepoint for forced start has not passed.
    if [ $FORCE_START -gt $LAST_TIME ] && [ $FORCE_START -le $CUR_TIME ]; then
      echo "[$(date +"%b %d %T")] Forced start triggered!"

      # Order missing batches according to amount of data.
      while read DIR; do
        nfastq=$(find $DIR/fastq -name *fastq | wc -l)
        nbases=$(( $nfastq * 4000 * 1050 ))

        echo -e $DIR"\t"$nbases
      done < /srv/rbd/covid19/processing/missing.txt |
      sort -r -k2,2 | 
      cut -f1 > tmp && mv tmp /srv/rbd/covid19/processing/missing.txt

      # Run workflow for all.
      while read DIR; do 
        if [ -z "$(pgrep -fl covid19_workflow.sh)" ]; then
          cd $(dirname $DIR) 
          #--------------------------------------------------------------------------------
          covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1 
          #mkdir $DIR/final_output
          #echo "hest" > $DIR/final_output/consensus.fasta
          #--------------------------------------------------------------------------------
          echo "[$(date +"%b %d %T")] workflow completed for $(basename $DIR)." 
        fi
      done < /srv/rbd/covid19/processing/missing.txt
    else

      ATTEMPT=0

      while read DIR; do 
    
        echo "[$(date +"%b %d %T")] Found a new batch: $(basename $DIR)"
  
        ### Test if there is enough data. If not, cancel the job.----------------------------------------------------------- 
        nfastq=$(find $DIR/fastq -name *.fastq | wc -l)
        nbases=$(( $nfastq * 4000 * 1050 ))
    
        # Grep the number Gb specified by user. If none, the default 20Gb are used.
        if [ -s $DIR/datacap.txt ]; then DATACAP=$(($(cat $DIR/datacap.txt) * 1000000000)); else DATACAP=20000000000; fi
    
        if [ $nbases -ge $DATACAP ]; then
         
          if [ -z "$(pgrep -fl covid19_workflow.sh)" ]; then
          
            # Run the workflow. 
            echo "[$(date +"%b %d %T")] $(basename $DIR) has enough data, continue with workflow processing." 
            cd $(dirname $DIR) 
            #--------------------------------------------------------------------------------
            covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1
            #mkdir $DIR/final_output
            #echo "hest" > $DIR/final_output/consensus.fasta
            #--------------------------------------------------------------------------------
            echo "[$(date +"%b %d %T")] workflow completed for $(basename $DIR)." 
          else 
            echo "[$(date +"%b %d %T")] covid19_workflow.sh is running elsewhere. Skipping and continue." 
          fi  
        fi
        sleep 3
      done < /srv/rbd/covid19/processing/missing.txt
    fi

    echo "[$(date +"%b %d %T")] Started waiting for new batches. Will check every 2 min." 

  else 
    ATTEMPT=$(($ATTEMPT+1))
  
    if [ ! $ATTEMPT -lt $LIFESIGN ]; then
      echo "[$(date +"%b %d %T")] Still waiting!"
      ATTEMPT=0
    fi
    sleep 2m
  fi

  LAST_TIME=$CUR_TIME

done