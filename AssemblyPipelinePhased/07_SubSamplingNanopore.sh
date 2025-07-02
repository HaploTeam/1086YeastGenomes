#!/bin/bash
#MSUB -r Phasing
#MSUB -n 71
#MSUB -c 1
#MSUB -T 86400
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

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
DOWNNA=DownSampleUnAssignedReads.cmd
MERGEFILES=MergeFiles.cmd
DOWNSAMPLE=DownSample.cmd
rm -f $DOWNNA
rm -f $MERGEFILES
rm -f $DOWNSAMPLE

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
WORKDIR=$BATCHDIR/07_DownsampleLR
NANOPOREDIR=$BATCHDIR/06_LRtagged

echo "--- Writting commands ---"

# Lower by 50% the coverage of unassigned reads
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	# Lower Coverage by 50%
	READS=$NANOPOREDIR/$S.NanoporeReads.haplotagged.noHPtag.fastq.gz
	OUT=$NANOPOREDIR/$S.NanoporeReads.haplotagged.noHPtag.Cov50percent.fastq.gz
	PERCENT=50
	echo "$FILTLONG --min_length 1000 --length_weight 10 --keep_percent $PERCENT --min_mean_q 9 $READS | gzip > $OUT" >> $DOWNNA
done

echo ""
echo "--- Running commands ---"
echo ""
ccc_mprun glost_launch $DOWNNA


# Merge unassigned reads with each haplotype
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		echo $S
		OUT=$WORKDIR/$S.$HP.TotalReads.fastq.gz
		echo "cat $NANOPOREDIR/$S.NanoporeReads.haplotagged.$HP.fastq.gz $NANOPOREDIR/$S.NanoporeReads.haplotagged.noHPtag.Cov50percent.fastq.gz > $OUT" >> $MERGEFILES
	done
done


echo ""
echo "--- Running commands ---"
echo ""
ccc_mprun glost_launch $MERGEFILES

# subsample at 30 and 40X
for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	for HP in HP1 HP2
	do
		# Compute coverage
		coverage=$(( $(zcat $WORKDIR/$S.$HP.TotalReads.fastq.gz | paste - - - - | cut -f 2 | tr -d '\n' | wc -c) / 12000000 ))
		if [ $coverage -ge 40 ]
		then
			COVERAGE=40
			OUT=$WORKDIR/$S.$HP.Cov$COVERAGE"X.SizePrior.minLength1kb.minQual9.fastq.gz"
			TARGET=$[ COVERAGE * 12000000 ]
			echo "$FILTLONG --min_length 1000 --length_weight 10 --target_bases $TARGET --min_mean_q 9 $WORKDIR/$S.$HP.TotalReads.fastq.gz | gzip > $OUT" >> $DOWNSAMPLE
		fi
		COVERAGE=30
		OUT=$WORKDIR/$S.$HP.Cov$COVERAGE"X.SizePrior.minLength1kb.minQual9.fastq.gz"
		TARGET=$[ COVERAGE * 12000000 ]
		echo "$FILTLONG --min_length 1000 --length_weight 10 --target_bases $TARGET --min_mean_q 9 $WORKDIR/$S.$HP.TotalReads.fastq.gz | gzip > $OUT" >> $DOWNSAMPLE
	done
done

echo ""
echo "--- Running commands ---"
echo ""
ccc_mprun glost_launch $DOWNSAMPLE

# Clean files
rm -f $MERGEFILES
rm -f $DOWNSAMPLE
rm -f $DOWNNA

