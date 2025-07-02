#!/bin/bash
#MSUB -r Polishing
#MSUB -n 128
#MSUB -c 8
#MSUB -Q long
#MSUB -T 259200
#MSUB -m scratch,work
#MSUB -q milan
#MSUB -o %N.%j.out 
#MSUB -e %N.%j.err

echo $(date)

module load extenv/fg
module load medaka
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/1_ASSEMBLY_Collapsed
DRAFTDIR=$BATCHDIR/10_RenamedAssemblies
WORKDIR=$BATCHDIR/12_Polishing

CMD="Polishing.cmd"
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch3.txt)
do
	DRAFT=$DRAFTDIR/$S.*.NoRedundantContigs.Nuclear.Renamed.fasta
	PREFIX=$(basename $DRAFT .fasta)

	NANOREADS=$BATCHDIR/01_SeqData/$S.NanoporeReads.fastq.gz

	# 1 round of Medaka
	echo "bash medaka.sh $DRAFT $NANOREADS $WORKDIR/$PREFIX.medakaPolished1.fasta" >> ${CMD}
done

ccc_mprun glost_launch ${CMD}
rm -f $CMD*

