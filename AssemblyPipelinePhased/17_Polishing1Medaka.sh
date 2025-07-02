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
module load medaka
module load glost

BATCHDIR=/ccc/scratch/cont007/fg0006/loeglerv/2_ASSEMBLY_Phased
DRAFTDIR=$BATCHDIR/16_RenamedAssemblies
WORKDIR=$BATCHDIR/18_Polishing

CMD="Polishing.cmd"
rm -f $CMD*

for S in $(cat $BATCHDIR/Batch3_Phased.txt)
do
	DRAFT=$DRAFTDIR/$S.PhasedGenome.fasta
	PREFIX=$(basename $DRAFT .fasta)

	NANOREADS=$BATCHDIR/01_SeqData/$S.NanoporeReads.fastq.gz

	# 1 round of Medaka
	echo "bash medaka.sh $DRAFT $NANOREADS $WORKDIR/$PREFIX.medakaPolished1.fasta" >> ${CMD}
done

ccc_mprun glost_launch ${CMD}
rm -f $CMD*

