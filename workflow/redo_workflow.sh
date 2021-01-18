#!/bin/bash

USAGE="$(basename "$0") [-h] [-i file -d flag -F flag] 

Arguments:
    -i  input text file, each line specifies full path to batch.
    -d  (flag) set to do dry-run, nothing is done. 
    -F  (flag) set to rerun complete workflow from demultiplexing onwards for all batches in -i, ignoring if workflow has been run on enough data. 

Usage:
It might be needed to rerun workflow for some batches, either because:
- The workflow was started before full data was reached
- The workflow needs to be rerun
This scripts controls if the workflow is rerun or not.
"
### Terminal Arguments ---------------------------------------------------------

# Import user arguments
while getopts ':hi:dF' OPTION; do
  case $OPTION in
    h) echo "$USAGE"; exit 1;;
    i) BATCHES=$OPTARG;;
    d) DRYRUN=true;;
    F) RERUN=true;;
    :) printf "missing argument for -$OPTARG\n" >&2; exit 1;;
    \?) printf "invalid option for -$OPTARG\n" >&2; exit 1;;
  esac
done

# Check missing arguments
MISSING="is missing but required. Exiting."
if [ -z ${BATCHES+x} ]; then echo "-i $MISSING"; exit 1; fi;

### CODE ---------------------------------------------------------------

# Reorder batches according to amount of data.
tfile=$(mktemp /tmp/redo.XXXXXXXXXX)

while read DIR; do
  nfastq=$(find $DIR/fastq -name *fastq | wc -l)
  nbases=$(( $nfastq * 4000 * 1050 ))

  echo -e $DIR"\t"$nbases

done < $BATCHES |
sort -rn -k2,2 > $tfile

# Takes as input the string 'rerun' to force-run workflow completely.
while read DIR GB; do 
 
  if [ "$RERUN" == "true" ]; then
    echo "[$(basename $DIR)] force runnning from scratch!"
    
    # Deleting output files. 
    if [ "$DRYRUN" != "true" ]; then
      rm -rf $DIR/demultiplexed
      rm -rf $DIR/processing
      rm -rf $DIR/QC
      rm -rf $DIR/final_output
    fi

    cd $(dirname $DIR) 
    echo "[$(basename $DIR)] Running workflow!"
    if [ "$DRYRUN" != "true" ]; then covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1; fi
  
  else
    if [ -s $DIR/data_used.txt ] && [ -s $DIR/datacap.txt ]; then
      CAP=$(cat $DIR/datacap.txt)
      USED=$(($(cat $DIR/data_used.txt) / 1000000000))
      AVAIL=$((GB / 1000000000))

      if [ $USED -ge $CAP ]; then
        echo "[$(basename $DIR)] This batch has been run with enough data, skipping."
      else 
        cd $(dirname $DIR) 
        echo "[$(basename $DIR)] This batch was run with ${USED}Gb data, ${AVAIL}Gb is now available. Running workflow!"
        if [ "$DRYRUN" != "true" ]; then covid19_workflow.sh -i $(basename $DIR) -s aau_long_v3.1; fi
      fi
    else
      echo "[$(basename $DIR)] Either data_used.txt and/or datacap.txt was not found, don't know what to do. You will have to manually handle this."
    fi
  fi
done < $tfile

rm $tfile