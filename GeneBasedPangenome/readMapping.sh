#!/bin/bash
# =============================================================================
# Project         : Full SaCe SNP matrix
# title           : 1_ReadsMapping.sh
# description     : This script will map illumina paired-end reads on a
#                   reference genome and output a BAM file. It is the first step
#                   for the construction of the Full SaCe SNP matrix. 
#                   The script requires 4 CPUs. 
# author          : vloegler
# date            : 2022/10/24
# version         : 2.0
# usage           : bash 1_ReadsMapping.sh $Strain $Reference $OutputPrefix $Reads_Forward $Reads_Reverse
# =============================================================================

module load bwa-mem2/2.2.1
module load samtools/1.15.1
module load gatk/4.2.3.0

STRAIN=$1 # Strain name
REF=$2 # Reference genome, Fasta format, indexed with bwa (bwa-mem2 index $REF)
PREFIX=$3 # Prefix of the outputted BAM file
BAM=$PREFIX.bam
ILLUMINAREADS1=$4 # Forward Illumina Reads, fastq or fastq.gz format (# Or single file containing reverse and forward)
ILLUMINAREADS2=$5 # Reverse Illumina Reads, fastq or fastq.gz format

# Map reads against reference genome
bwa-mem2 mem -t 20 -U 0 -L 0,0 -O 4,4 -T 20 $REF $ILLUMINAREADS1 $ILLUMINAREADS2 | samtools sort -o $BAM -T reads.$(basename $PREFIX).tmp
# Add read groups
gatk AddOrReplaceReadGroups -I $BAM -O $BAM.ReadGroups --RGID $STRAIN --RGLB $STRAIN --RGPL ILLUMINA --RGPU $STRAIN --RGSM $STRAIN
rm -f $BAM
mv $BAM.ReadGroups $BAM
# Index BAM
samtools index $BAM
