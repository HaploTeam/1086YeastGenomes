#!/bin/bash
# This script run 1 round of Hapo-G short reads polishing
# Input arguments:
# hapoG.sh $DraftGenomeAssembly.fasta $ShortReads_1.fastq.gz $ShortReads_2.fastq.gz $OutputFile

echo $(date)

module load extenv/fg
module load python
module load samtools
module load bwa
module load hapog

WORKDIR=$(pwd)
RAWASSEMBLY=$1
SHORTREADS1=$2
SHORTREADS2=$3
SHORTREADS="$SHORTREADS1 $SHORTREADS2"
OUTPUT=$4

# Run Hapo-G
hapog --genome $RAWASSEMBLY --pe1 $SHORTREADS1 --pe2 $SHORTREADS2 -o $WORKDIR/$(basename $RAWASSEMBLY)_HapoG -t 24

# Clean output files
cp $WORKDIR/$(basename $RAWASSEMBLY)_HapoG/hapog_results/hapog.fasta $OUTPUT
rm -rf $WORKDIR/$(basename $RAWASSEMBLY)_HapoG

