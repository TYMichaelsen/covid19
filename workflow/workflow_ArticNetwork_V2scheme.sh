#/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log.out 2>&1
cd /srv/rbd
singularity shell docker://kasperskytte/covid19
cd rhk/test2/

#eval "$(conda shell.bash hook)"
source activate artic-ncov2019-medaka
# Settings
FASTQ=/srv/rbd/covid19/data/COVID19/analysis/AAU-batch001/data/aaubatch001.fastq
RUNID=$(basename $FASTQ | sed 's/.fastq//')
SCHEMEDIR=/srv/rbd/covid19/data/COVID19/artic-ncov2019/primer_schemes
THREADS=120
REF=/srv/rbd/covid19/data/COVID19/artic-ncov2019/primer_schemes/nCoV-2019/V2/nCoV-2019.reference.fasta

# setup output folders
mkdir -p data
mkdir -p tempo/TMPDIR/
mkdir -p tempo/trim
mkdir -p tempo/demultiplexed/
mkdir -p tempo/articminion/
mkdir -p results/genomes/
mkdir -p results/coverage/
mkdir -p results/mapped_fastq/
mkdir -p results/N_counts/

#################################################
# Basecall reads	 							#
#################################################

#################################################
# Pre-process samples 							#
#################################################
# Filter fastq by size (https://www.biostars.org/p/66996/)
OUTPUTFILE=tempo/filtered.fastq
if [ -s $OUTPUTFILE ]; then rm $OUTPUTFILE; awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 400 && length(seq) <= 700) {print header, seq, qheader, qseq}}' < $FASTQ > tempo/filtered.fastq; 
else
awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 400 && length(seq) <= 700) {print header, seq, qheader, qseq}}' < $FASTQ > tempo/filtered.fastq
fi
# Run porechop in parallel
TRIM_DIR=tempo/trim
rm -r $TRIM_DIR
mkdir -p $TRIM_DIR
cat tempo/filtered.fastq | parallel --progress -j $THREADS -L 4 --round-robin --pipe --tmpdir tempo/TMPDIR \
    "cat > $TRIM_DIR/{#}.tmp;\
    porechop \
      -i $TRIM_DIR/{#}.tmp \
      --threads 1 \
	  --verbosity 2 \
	  --untrimmed \
	  --barcode_dir $TRIM_DIR/{#} \
	  --native_barcodes \
	  --discard_middle \
	  --require_two_barcodes \
	  --barcode_threshold 80 \
	  --check_reads 10000 >> out 2>>error"

# Concatanate files
find tempo/trim/ -name "*.fastq" | while read file
do
        base=${file##*/}    # ditch the directory portion of the filename
        cat $file >> tempo/demultiplexed/$RUNID"_"$base
done

rm tempo/demultiplexed/$RUNID"_"none.fastq

#################################################
# Process samples 								#
#################################################
# Run reference "assemblies"
FILES=tempo/demultiplexed/*.fastq
for f in $FILES
do
SAMPLE=$(basename $f | sed 's/.fastq//')
OUTPUTFILE=results/coverage/$SAMPLE.cov.tsv
if [ -s $OUTPUTFILE ]; then echo "$OUTPUTFILE has already been generated"; 
else
# Basic artic minion workflow
artic minion --medaka --minimap2 --normalise 200 --threads $THREADS --scheme-directory $SCHEMEDIR --read-file tempo/demultiplexed/$SAMPLE.fastq nCoV-2019/V2 $SAMPLE
mv $SAMPLE* tempo/articminion/
# Generate coverage
echo -e 'scaffold\tposition\tcoverage' > results/coverage/$SAMPLE.cov.tsv
samtools depth -a tempo/articminion/$SAMPLE.sorted.bam >> results/coverage/$SAMPLE.cov.tsv

mv tempo/articminion/$SAMPLE.consensus.fasta results/genomes/
# Extract mapped reads
samtools fastq --threads $THREADS -F 4 tempo/articminion/$SAMPLE.sorted.bam > results/mapped_fastq/$SAMPLE"_virus".fastq
# Calculate the number of Ns in the sequences
awk '!/^>/ { next } { getline seq; len=length(seq); Nn=gsub(/N/,"",seq) } {sub(/^>/,"",$0); print $0 "\t" Nn}' results/genomes/$SAMPLE.consensus.fasta > results/N_counts/$SAMPLE"N_count.tsv"

fi
done
conda deactivate

#################################################
# Plotting		 								#
#################################################
