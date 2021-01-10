#!/bin/bash

BATCHES=$1
RERUN=$2

# It might be needed to rerun workflow for some batches, either because:
# - The workflow was started before full data was reached
# - The workflow needs to be rerun

# Reorder batches according to amount of data.
BATCHES_ORD=$(while read DIR; do
  nfastq=$(find $DIR/fastq -name *fastq | wc -l)
  nbases=$(( $nfastq * 4000 * 1050 ))

  echo -e $DIR"\t"$nbases

done < $BATCHES |
sort -rn -k2,2 | 
cut -f1 -)

# Takes as input the string 'rerun' to force-run workflow completely.
for DIR in $BATCHES_ORD; do 
 
  if [ ! -z "$RERUN" ]; then
    echo "[$(basename $DIR)] force runnning from scratch!"
     
    rm -rf $DIR/demultiplexed
    rm -rf $DIR/processing
    rm -rf $DIR/QC
    rm -rf $DIR/final_output

    cd $(dirname $DIR) 
    echo "[$(basename $DIR)] Running workflow!"
    covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1

  else
    if [ -s $DIR/data_used.txt ] && [ -s $DIR/datacap.txt ]; then
      CAP=$(($(cat $DIR/datacap.txt) * 1000000000))
      USED=$(cat $DIR/data_used.txt)
      if [ $USED -ge $CAP ]; then
        echo "[$(basename $DIR)] This batch has been run with enough data, skipping."
      else 
        cd $(dirname $DIR) 
        echo "[$(basename $DIR)] Running workflow!"
	covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1
      fi
    else
      echo "[$(basename $DIR)] Either data_used.txt and/or datacap.txt was not found, don't know what to do. You will have to manually handle this."
    fi
  fi
done 
