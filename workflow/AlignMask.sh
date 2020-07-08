#!/usr/bin/env bash
OUTDIR=$1
SEQS=$2
REF=$3
DEPEND_DIR=$4
THREADS=$5

source activate nextstrain

# Do alignment in chunks.
rm -rf $OUTDIR/split_fasta; mkdir $OUTDIR/split_fasta
rm -rf $OUTDIR/split_align; mkdir $OUTDIR/split_align

split -d -l 40 --additional-suffix=.fasta $OUTDIR/$SEQS $OUTDIR/split_fasta/
FASTAFILES=`ls $OUTDIR/split_fasta` 

parallel -j$THREADS \
         '
augur align \
--sequences {2}/split_fasta/{1} \
--reference-sequence {3} \
--output {2}/split_align/{1}.aligned \
--nthreads 1 \
--remove-reference \
--fill-gaps &> {2}/log.out

' ::: $FASTAFILES ::: $OUTDIR ::: $REF

cat $OUTDIR/split_align/*.aligned > $OUTDIR/aligned.fasta
cat $OUTDIR/split_align/*.log > $OUTDIR/aligned.log

rm -r $OUTDIR/split_fasta
rm -r $OUTDIR/split_align

# Mask bases.
mask_sites="13402 24389 24390"

python3 ${DEPEND_DIR}/mask-alignment.py \
        --alignment $OUTDIR/aligned.fasta \
        --mask-from-beginning 100 \
        --mask-from-end 50 \
        --mask-sites $mask_sites \
        --output $OUTDIR/masked.fasta