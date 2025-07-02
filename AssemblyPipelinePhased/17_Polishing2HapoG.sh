#!/bin/bash
#MSUB -r Polishing
#MSUB -n 64
#MSUB -c 8
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load extenv/fg
module load python
module load r
module load nanoporetech
module load samtools
module load racon
module load minimap2
module load samtools
module load bwa
module load hapog
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
DRAFTDIR=$BATCHDIR/18_Polishing
WORKDIR=$BATCHDIR/18_Polishing

CMD="Polishing.cmd"
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch2_Phased.txt)
do
	DRAFT=$DRAFTDIR/$S.*.medakaPolished1.fasta
	PREFIX=$(basename $DRAFT .fasta)
	ILLUREADS="$BATCHDIR/17_IlluminaSeqData/$S.illuminaSeqData_1.fastq.gz $BATCHDIR/17_IlluminaSeqData/$S.illuminaSeqData_2.fastq.gz"

	# 1 round of hapoG
	echo "bash hapoG.sh $DRAFT $ILLUREADS $WORKDIR/$PREFIX.hapoGPolished1.fasta" >> ${CMD}6

done
ccc_mprun glost_launch ${CMD}6
rm -f $CMD*

