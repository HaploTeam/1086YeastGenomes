#!/bin/bash
# This script run 1 round of Racon long reads polishing
# Input arguments:
# racon.sh $DraftGenomeAssembly.fasta $LongReads.fastq.gz $OutputFile.fasta

module load extenv/fg
module load python
module load r
module load samtools
module load racon
module load minimap2

WORKDIR=$(pwd)
RAWASSEMBLY=$1
NANOPOREREADS=$2
OUTPUT=$3
MAPPINGDIR=$WORKDIR/$(basename $RAWASSEMBLY)_NanoporeMapping
mkdir $MAPPINGDIR
SAM=$MAPPINGDIR/$(basename $RAWASSEMBLY)

# Map Nanopore reads to the draft assembly
minimap2 -ax map-ont -t 20 $RAWASSEMBLY $NANOPOREREADS | samtools sort -o $SAM.sam -T reads.$(basename $RAWASSEMBLY).tmp

# Polish with Racon
racon -t 20 $NANOPOREREADS $SAM.sam $RAWASSEMBLY > $OUTPUT

# Clean files
rm -rf $MAPPINGDIR

