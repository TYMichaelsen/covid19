#!/bin/bash
# By Thomas Y. Michaelsen

cd /srv/rbd/test_workflow

mkdir -p analysis/data 

# Sweep the export folder for genomes.
cat ../export/*_export/sequences.fasta > analysis/data/sequences.fasta

# Export the associated metadata.
