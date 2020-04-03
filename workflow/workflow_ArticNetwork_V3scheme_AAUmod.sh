
#/bin/bash
exec 3>&1 4>&2
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1>log.out 2>&1
cd /srv/rbd
singularity shell docker://kasperskytte/covid19
cd rhk/test3/

#eval "$(conda shell.bash hook)"
source activate artic-ncov2019-medaka
# Settings
FASTQDIR=/srv/rbd/covid19/data/CJ023/demultiplexed/
RUNID=CJ023
/srv/rbd/covid19/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/
SCHEMEDIR=/srv/rbd/covid19/software/artic-ncov2019/primer_schemes/
THREADS=120
REF=/srv/rbd/covid19/software/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta

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

cp $FASTQDIR/*.fastq tempo/demultiplexed/

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
artic minion --medaka --minimap2 --normalise 200 --threads $THREADS --scheme-directory $SCHEMEDIR --read-file tempo/demultiplexed/$SAMPLE.fastq nCoV-2019/V3 $SAMPLE
mv $SAMPLE* tempo/articminion/
# Generate coverage
echo -e 'scaffold\tposition\tcoverage' > results/coverage/$SAMPLE.cov.tsv
samtools depth -a tempo/articminion/$SAMPLE.sorted.bam >> results/coverage/$SAMPLE.cov.tsv

mv tempo/articminion/$SAMPLE.consensus.fasta results/genomes/
# Extract mapped reads
samtools fastq --threads $THREADS -F 4 tempo/articminion/$SAMPLE.sorted.bam > results/mapped_fastq/$SAMPLE"_virus".fastq
md5sum results/mapped_fastq/$SAMPLE"_virus".fastq >> results/md5sums.txt
# Calculate the number of Ns in the sequences
awk '!/^>/ { next } { getline seq; len=length(seq); Nn=gsub(/N/,"",seq) } {sub(/^>/,"",$0); print $0 "\t" Nn}' results/genomes/$SAMPLE.consensus.fasta > results/N_counts/$SAMPLE"N_count.tsv"

fi
done
conda deactivate

#################################################
# Plotting		 								#
#################################################
