#!/bin/bash
#MSUB -r Filtlong
#MSUB -n 256
#MSUB -c 1
#MSUB -T 86400
#MSUB -q milan
#MSUB -m scratch,work
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load extenv/fg
module load python
module load r
module load glost

# =======================================================================================
# This script downsample Nanopore reads without using an illumina reference with FiltLong
# =======================================================================================

# Software path
FILTLONG=/ccc/work/cont007/fg0006/loeglerv/Soft/Filtlong-0.2.1/bin/filtlong
# Commands to run in parallel
DOWNSAMPLE=DownSample.cmd
rm -f $DOWNSAMPLE

NANOPOREDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed/01_SeqData

echo "--- Writting commands ---"



# Subsample at 40X
for S in $(cat Batch2_Down40X.txt)
do
	COVERAGE=40
	IN=$NANOPOREDIR/$S.NanoporeReads.fastq.gz
	OUT=$NANOPOREDIR/$S.Cov$COVERAGE"X.SizePrior.minLength1kb.minQual9.fastq.gz"
	TARGET=$[ COVERAGE * 12000000 ]
	echo "$FILTLONG --min_length 1000 --length_weight 10 --target_bases $TARGET --min_mean_q 9 $IN | gzip > $OUT" >> $DOWNSAMPLE
done

# Subsample at 30X
for S in $(cat Batch2.txt)
do
	COVERAGE=30
	IN=$NANOPOREDIR/$S.NanoporeReads.fastq.gz
	OUT=$NANOPOREDIR/$S.Cov$COVERAGE"X.SizePrior.minLength1kb.minQual9.fastq.gz"
	TARGET=$[ COVERAGE * 12000000 ]
	echo "$FILTLONG --min_length 1000 --length_weight 10 --target_bases $TARGET --min_mean_q 9 $IN | gzip > $OUT" >> $DOWNSAMPLE
done

echo ""
echo "--- Running commands ---"
echo ""
ccc_mprun glost_launch $DOWNSAMPLE

# Clean files
rm -f $DOWNSAMPLE


