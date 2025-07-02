#!/bin/bash
# This script run 1 round of Medaka long reads polishing
# Input arguments:
# medaka.sh $DraftGenomeAssembly.fasta $LongReads.fastq.gz $OutputFile.fasta

echo $(date)

module load extenv/fg
module load medaka

WORKDIR=$(pwd)
RAWASSEMBLY=$1
NANOPOREREADS=$2
OUTPUT=$3

NPROC=20
BASECALLS=$NANOPOREREADS
RUN=r941_prom_sup_g507

medaka_consensus -i $BASECALLS -d $RAWASSEMBLY -o $WORKDIR/$(basename $RAWASSEMBLY)_DIR -t $NPROC -m $RUN
mv $WORKDIR/$(basename $RAWASSEMBLY)_DIR/consensus.fasta $OUTPUT
rm -rf $WORKDIR/$(basename $RAWASSEMBLY)_DIR
rm -f $RAWASSEMBLY.fai
rm -f $RAWASSEMBLY.mmi

